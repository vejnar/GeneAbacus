//
// Copyright © 2015 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package feature

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"hash/adler32"
	"math"
	"os"
	"strconv"
	"strings"

	"github.com/pierrec/lz4"

	"git.sr.ht/~vejnar/GeneAbacus/lib/cmapper"
)

const (
	bedGraphPrecision = 0.000001
)

type FeatureExt struct {
	*Feature
	CoordMapper *cmapper.CoordMapper
	Counts      []float64
	Profile     []float32
}

func ExtendFeatures(features []Feature, countMultis []int, doProfile bool, profileOverhang int) ([]*FeatureExt, error) {
	featureExts := make([]*FeatureExt, len(features))
	for ifeat := 0; ifeat < len(features); ifeat++ {
		// New
		fe := FeatureExt{Feature: &features[ifeat]}
		if fe.ID != uint32(ifeat) {
			return featureExts, fmt.Errorf("Wrong feature ID")
		}
		// Init. count
		fe.Counts = make([]float64, 1+len(countMultis)*2)
		// Length
		fe.Counts[0] = float64(IntervalsLength(fe.Coords))
		// Init. profile
		if doProfile {
			// Deep-copy
			coords := make([][]int, len(fe.Coords))
			for i := 0; i < len(fe.Coords); i++ {
				coords[i] = make([]int, len(fe.Coords[i]))
				copy(coords[i], fe.Coords[i])
			}
			// Add overhang
			coords[0][0] -= profileOverhang
			coords[len(coords)-1][1] += profileOverhang
			// CoordMapper
			fe.CoordMapper = &cmapper.CoordMapper{CoordsGenome: coords, Strand: fe.Strand}
			fe.CoordMapper.Init()
			// Profile
			fe.Profile = make([]float32, fe.CoordMapper.Length)
		}
		// Append feature
		featureExts[ifeat] = &fe
	}
	return featureExts, nil
}

func WriteCounts(featureExts []*FeatureExt, countPath string, countMultis []int, totals []float64, appendOutput bool) error {
	// Append or Create flag
	var fg int
	if appendOutput {
		fg = os.O_APPEND | os.O_CREATE | os.O_WRONLY
	} else {
		fg = os.O_RDWR | os.O_CREATE | os.O_TRUNC
	}
	if f, err := os.OpenFile(countPath, fg, 0666); err != nil {
		return err
	} else {
		// Number of commas to add
		ncomma := len(countMultis) - 1
		// Write header
		f.WriteString("\"name\",\"length\",")
		for i, cm := range countMultis {
			fmt.Fprintf(f, "\"count_%d\",\"rpkm_%d\"", cm, cm)
			if i < ncomma {
				f.WriteString(",")
			}
		}
		f.WriteString("\n")
		// Totals
		ncomma = len(countMultis) * 2
		f.WriteString("\"total\",")
		for i, t := range totals {
			f.WriteString(strconv.FormatFloat(t, 'f', -1, 64))
			if i < ncomma {
				f.WriteString(",")
			}
		}
		f.WriteString("\n")
		// Write counts
		for _, feat := range featureExts {
			fmt.Fprintf(f, "\"%s\",", feat.Name)
			for i, c := range feat.Counts {
				f.WriteString(strconv.FormatFloat(c, 'f', -1, 32))
				if i < ncomma {
					f.WriteString(",")
				}
			}
			f.WriteString("\n")
		}
		// Close CSV
		f.Close()
	}
	return nil
}

type GenericWriter interface {
	Write(buf []byte) (n int, err error)
	Close() error
}

func WriteProfiles(featureExts []*FeatureExt, featuresMapping map[string]string, profilePath string, profileFormat string, appendOutput bool) error {
	var profileZip string
	var mapName bool
	if len(featuresMapping) > 0 {
		mapName = true
	}
	if strings.Contains(profileFormat, "+") {
		doubleFormat := strings.Split(profileFormat, "+")
		profileFormat, profileZip = doubleFormat[0], doubleFormat[1]
	}
	// Append or Create flag
	var fg int
	if appendOutput {
		fg = os.O_APPEND | os.O_CREATE | os.O_WRONLY
	} else {
		fg = os.O_RDWR | os.O_CREATE | os.O_TRUNC
	}
	if f, err := os.OpenFile(profilePath, fg, 0666); err != nil {
		return err
	} else {
		var writer GenericWriter
		switch profileZip {
		case "lz4":
			writer = lz4.NewWriter(f)
		case "lz4hc":
			lzWriter := lz4.NewWriter(f)
			lzWriter.Header = lz4.Header{CompressionLevel: 9}
			writer = lzWriter
		default:
			writer = f
		}
		switch profileFormat {
		case "bedgraph":
			var diff float64
			for _, feat := range featureExts {
				var stepStart int
				var stepValue, currentValue float32
				var name string
				if mapName {
					name = MapName(feat.Name, featuresMapping)
				} else {
					name = feat.Name
				}
				for ip := 0; ip < len(feat.Profile); ip++ {
					currentValue = feat.Profile[ip]
					if diff = math.Abs(float64(currentValue - stepValue)); diff > bedGraphPrecision {
						if stepValue != 0. {
							fmt.Fprintf(writer, "%s\t%d\t%d\t%f\n", name, stepStart, ip, stepValue)
						}
						stepStart = ip
						stepValue = currentValue
					}
				}
			}
		case "binary":
			// Version
			var version uint8
			version = 3
			binary.Write(writer, binary.LittleEndian, version)
			// Features and total lengths
			var totalLength, l uint32
			bufChecksum := new(bytes.Buffer)
			for i := 0; i < len(featureExts); i++ {
				l = uint32(featureExts[i].Length())
				err := binary.Write(bufChecksum, binary.LittleEndian, l)
				if err != nil {
					return err
				}
				totalLength += l
			}
			// Write total length
			binary.Write(writer, binary.LittleEndian, totalLength)
			// Checksum
			checksum := adler32.Checksum(bufChecksum.Bytes())
			err = binary.Write(writer, binary.LittleEndian, checksum)
			if err != nil {
				return err
			}
			// Write profiles
			for _, feat := range featureExts {
				err = binary.Write(writer, binary.LittleEndian, feat.Profile)
				if err != nil {
					return err
				}
			}
		case "csv":
			for _, feat := range featureExts {
				var name string
				if mapName {
					name = MapName(feat.Name, featuresMapping)
				} else {
					name = feat.Name
				}
				fprofile := fmt.Sprintf("%v", feat.Profile)
				fmt.Fprintf(writer, "%s,%d,%s\n", name, len(feat.Profile), fprofile[1:len(fprofile)-1])
			}
		}
		writer.Close()
		f.Close()
	}
	return nil
}
