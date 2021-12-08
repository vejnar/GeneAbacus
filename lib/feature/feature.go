//
// Copyright (C) 2015-2021 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package feature

import (
	"bufio"
	"encoding/json"
	"fmt"
	"os"
	"strconv"
	"strings"
)

type Feature struct {
	ID     uint32
	Name   string
	Chrom  string
	Strand int8
	Coords [][]int
}

// Length returns the length of feature
func (feat Feature) Length() (length int) {
	for _, coords := range feat.Coords {
		length += coords[1] - coords[0]
	}
	return
}

// Sorting functions: By Name
// Use it with: sort.Sort(feature.ByName(features))
type ByName []Feature

func (f ByName) Len() int           { return len(f) }
func (f ByName) Swap(i, j int)      { f[i], f[j] = f[j], f[i] }
func (f ByName) Less(i, j int) bool { return f[i].Name < f[j].Name }

// Sorting functions: By Chrom
type ByChrom []Feature

func (f ByChrom) Len() int           { return len(f) }
func (f ByChrom) Swap(i, j int)      { f[i], f[j] = f[j], f[i] }
func (f ByChrom) Less(i, j int) bool { return f[i].Chrom < f[j].Chrom }

// OpenFON parses a "Feature Object Notation" string and returns a list of Feature
func OpenFON(jpath, fonName, fonChrom, fonStrand, fonCoords string) (features []Feature, err error) {
	jfos, err := os.Open(jpath)
	if err != nil {
		return
	}
	defer jfos.Close()

	d := json.NewDecoder(jfos)
	d.UseNumber()
	var rawFON interface{}
	if err = d.Decode(&rawFON); err != nil {
		err = fmt.Errorf("Error while parsing JSON feature file %s", jpath)
		return
	}

	// FON
	fon := rawFON.(map[string]interface{})

	// FON version
	var version int64
	if version, err = fon["fon_version"].(json.Number).Int64(); err != nil {
		return
	} else if version != 1 {
		err = fmt.Errorf("Unknown FON version %d", version)
		return
	}

	// Get features
	rawFeatures := fon["features"].([]interface{})
	for i, rf := range rawFeatures {
		mf := rf.(map[string]interface{})
		// Strand
		var istrand int8
		strand := mf[fonStrand].(string)
		if strand == "+" {
			istrand = 1
		} else if strand == "-" {
			istrand = -1
		}
		f := Feature{ID: uint32(i), Name: mf[fonName].(string), Chrom: mf[fonChrom].(string), Strand: istrand}
		// Add coordinates
		f.Coords = make([][]int, len(mf[fonCoords].([]interface{})))
		for j, cj := range mf[fonCoords].([]interface{}) {
			f.Coords[j] = make([]int, 2)
			for k, ck := range cj.([]interface{}) {
				n, _ := ck.(json.Number).Int64()
				f.Coords[j][k] = int(n)
			}
		}
		features = append(features, f)
	}
	return
}

// OpenTAB parses a two column tabulated files with name and length of feature and returns a list of Feature
func OpenTAB(tpath string, strand int8) (features []Feature, err error) {
	tfos, err := os.Open(tpath)
	if err != nil {
		return
	}
	defer tfos.Close()

	var i uint32
	var length int
	tscanner := bufio.NewScanner(tfos)
	for tscanner.Scan() {
		fields := strings.Split(tscanner.Text(), "\t")
		length, err = strconv.Atoi(fields[1])
		if err != nil {
			return
		}
		f := Feature{ID: i, Name: fields[0], Chrom: fields[0], Strand: strand, Coords: [][]int{[]int{0, length}}}
		features = append(features, f)
		i++
	}
	if err = tscanner.Err(); err != nil {
		return
	}
	return
}

// IntervalsLength returns the length covered by all intervals (0-based [start,end))
func IntervalsLength(intervals [][]int) (length int) {
	for _, iv := range intervals {
		length += iv[1] - iv[0]
	}
	return
}
