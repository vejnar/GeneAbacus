//
// Copyright (C) 2015-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package main

import (
	"context"
	"fmt"
	"io"
	"math/rand"
	"os"
	"os/exec"
	"strconv"
	"time"

	"git.sr.ht/~vejnar/GeneAbacus/lib/esam"
	"git.sr.ht/~vejnar/GeneAbacus/lib/feature"
	"git.sr.ht/~vejnar/GeneAbacus/lib/profile"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/biogo/store/interval"

	"golang.org/x/sync/errgroup"

	"gopkg.in/fatih/set.v0"
)

const (
	cacheLength        = 2
	cacheProfileLength = 100
	sPairLength        = 10
)

type Packet struct {
	ID             uint32
	Counts         []float64
	ProfileChanges *profile.ProfileChange
}

type Cache struct {
	Packets     []Packet
	LastPacket  int
	InputCount  float64
	MultiCounts []float64
}

func NewCache(size int, nMulti int) *Cache {
	c := Cache{}
	c.MultiCounts = make([]float64, nMulti)
	c.Packets = make([]Packet, size)
	for i := 0; i < size; i++ {
		// Count
		c.Packets[i].Counts = make([]float64, nMulti)
		// Profile
		c.Packets[i].ProfileChanges = profile.NewProfileChange(cacheProfileLength)
	}
	return &c
}

func (c *Cache) Grow() {
	osize := len(c.Packets)
	nsize := int(float64(osize) * 1.5)
	c.Packets = append(c.Packets, make([]Packet, nsize-osize)...)
	for i := osize; i < nsize; i++ {
		// Count
		c.Packets[i].Counts = make([]float64, len(c.Packets[0].Counts))
		// Profile
		c.Packets[i].ProfileChanges = profile.NewProfileChange(cacheProfileLength)
	}
}

type Pair struct {
	Reads     []*sam.Record
	OnlyRead1 bool
}

// AddCommas adds commas after every 3 characters.
func AddCommas(s string) string {
	if len(s) <= 3 {
		return s
	} else {
		return AddCommas(s[0:len(s)-3]) + "," + s[len(s)-3:]
	}
}

func Max(x, y int) int {
	if x > y {
		return x
	}
	return y
}

func Abs(n int) int {
	if n < 0 {
		return -n
	}
	return n
}

func OpenSAM(pathSAM esam.PathSAM, cmd []string, nWorker1 int) (f *os.File, pp io.ReadCloser, rr sam.RecordReader, err error) {
	if pathSAM.Binary {
		f, err = os.Open(pathSAM.Path)
		if err != nil {
			return f, pp, rr, err
		}
		rr, err = bam.NewReader(f, nWorker1)
	} else {
		if len(cmd) == 0 {
			f, err = os.Open(pathSAM.Path)
			if err != nil {
				return f, pp, rr, err
			}
			rr, err = sam.NewReader(f)
		} else {
			cmd = append(cmd, pathSAM.Path)
			p := exec.Command(cmd[0], cmd[1:]...)
			if pp, err = p.StdoutPipe(); err != nil {
				return f, pp, rr, err
			}
			if err = p.Start(); err != nil {
				return f, pp, rr, err
			}
			rr, err = sam.NewReader(pp)
		}
	}
	return f, pp, rr, nil
}

func GetSAMHeader(pathSAM esam.PathSAM, cmd []string) (*sam.Header, error) {
	var err error
	f, err := os.Open(pathSAM.Path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	if pathSAM.Binary {
		rr, err := bam.NewReader(f, 1)
		if err != nil {
			return nil, err
		}
		return rr.Header(), nil
	} else {
		if len(cmd) == 0 {
			f, err := os.Open(pathSAM.Path)
			if err != nil {
				return nil, err
			}
			rr, err := sam.NewReader(f)
			return rr.Header(), nil
		} else {
			var pp io.ReadCloser
			var err error
			cmd = append(cmd, pathSAM.Path)
			p := exec.Command(cmd[0], cmd[1:]...)
			if pp, err = p.StdoutPipe(); err != nil {
				return nil, err
			}
			defer pp.Close()
			if err := p.Start(); err != nil {
				return nil, err
			}
			rr, err := sam.NewReader(pp)
			return rr.Header(), nil
		}
	}
	return nil, nil
}

func PConFeature(pathSAMs []esam.PathSAM, SAMCmdIn []string, features []feature.Feature, featuresMapping map[string]string, trees map[string]map[int8]*interval.IntTree, readLengths []int, fragmentMinLength int, fragmentMaxLength int, randProportion float32, paired bool, libraryR1Strand int8, ignoreNHTag bool, inProperPair bool, minMappingQuality byte, minOverlap int, countMultis []int, countTotals []float64, countTotalInput bool, countTotalRealRead bool, countInProfile bool, countPath string, profileType int, profileMulti int, profileOverhang int, profileNoCoordMapping bool, profileUntemplated int, profileNoUntemplated bool, profileExtensionLength int, profilePositionFraction float64, profileNorm bool, profileMultiTotalCol int, profilePaths []string, profileFormats []string, appendOutput bool, pathReport string, pathSAMOut esam.PathSAM, nWorker int, timeStart time.Time, verboseLevel int) (nAlign uint64, err error) {
	// Compute profile(s) ?
	var doProfile bool
	if profileType != profile.ProfileTypeNone {
		doProfile = true
	}
	// Workers
	nWorker1 := Max(1, int(nWorker/2.))
	nWorker2 := Max(1, nWorker-nWorker1)

	// Init. extended features
	var featureExts []*feature.FeatureExt
	featureExts, err = feature.ExtendFeatures(features, countMultis, doProfile, profileOverhang)
	if err != nil {
		return nAlign, err
	}

	// Init. input counter
	var inputCount float64

	// Init. read or multiplicity counter
	var multiSets []set.Interface
	var multisCounts []float64
	if countTotalRealRead {
		if verboseLevel > 0 {
			timeNow := time.Now()
			fmt.Printf("%.1fmin - Init. total full-read counter\n", timeNow.Sub(timeStart).Minutes())
		}
		for icm := 0; icm < len(countMultis); icm++ {
			multiSets = append(multiSets, set.New(set.ThreadSafe))
		}
	} else {
		if verboseLevel > 0 {
			timeNow := time.Now()
			fmt.Printf("%.1fmin - Init. total proportion-read counter\n", timeNow.Sub(timeStart).Minutes())
		}
		multisCounts = make([]float64, len(countMultis))
	}

	// Open output SAM
	var doOutSAM bool
	var samWriter *sam.Writer
	if pathSAMOut.Path != "" {
		// Open output file
		f, err := os.Create(pathSAMOut.Path)
		if err != nil {
			return nAlign, err
		}
		// Get SAM header
		samHeader, err := GetSAMHeader(pathSAMs[0], SAMCmdIn)
		if err != nil {
			return nAlign, err
		}
		// Create SAM writer
		samWriter, err = sam.NewWriter(f, samHeader, sam.FlagDecimal)
		if err != nil {
			return nAlign, err
		}
		defer f.Close()
		doOutSAM = true
	}

	// Init context
	ctx, _ := context.WithCancel(context.Background()) // cancelFunc not used
	// Start sync errgroup
	g, gctx := errgroup.WithContext(ctx)

	// Start receiving channel
	chFinal := make(chan *Cache, nWorker*10)
	// Start alignment channel
	chAln := make(chan []*Pair, nWorker*10)

	//go func() {
	g.Go(func() error {
		defer close(chAln)
		timeLog := time.Now()
		for _, pathSAM := range pathSAMs {
			var f *os.File
			var pp io.ReadCloser
			var rr sam.RecordReader
			var err error
			var iPair int
			sPair := make([]*Pair, sPairLength)
			if verboseLevel > 0 {
				timeNow := time.Now()
				fmt.Printf("%.1fmin - Opening %s\n", timeNow.Sub(timeStart).Minutes(), pathSAM.Path)
			}
			// Open SAM
			f, pp, rr, err = OpenSAM(pathSAM, SAMCmdIn, nWorker1)
			if err != nil {
				return err
			}
			defer f.Close()
			if pp != nil {
				defer pp.Close()
			}

			// Loop over reads
			var isRead1First, isRead2First, read1Mapped, mateMapped bool
			for {
				// Next read
				var pair Pair
				var aread, areadM *sam.Record
				aread, err = rr.Read()
				if err == io.EOF {
					break
				} else if err != nil {
					return err
				}
				read1Mapped = aread.Flags&sam.Unmapped == 0
				// Get mate
				if paired {
					mateMapped = aread.Flags&sam.MateUnmapped == 0
					if mateMapped {
						for {
							areadM, err = rr.Read()
							if err == io.EOF {
								break
							} else if err != nil {
								return err
							}
							// If alignment is not supplementary, let's keep it
							if areadM.Flags&sam.Supplementary == 0 {
								break
							}
						}
						// Check read1 and read2 names are the same
						if aread.Name != areadM.Name {
							return fmt.Errorf("Differerent names for Read1 %s and Read2 %s", aread.Name, areadM.Name)
						}
					}
					// Combine read(s) into pair
					isRead1First = aread.Flags&sam.Read1 != 0
					isRead2First = aread.Flags&sam.Read2 != 0
					if read1Mapped && mateMapped {
						if isRead1First {
							pair.Reads = append(pair.Reads, aread, areadM)
						} else {
							pair.Reads = append(pair.Reads, areadM, aread)
						}
					} else if read1Mapped {
						pair.Reads = append(pair.Reads, aread)
						if isRead1First {
							pair.OnlyRead1 = true
						}
					} else if mateMapped {
						pair.Reads = append(pair.Reads, areadM)
						if isRead2First {
							pair.OnlyRead1 = false
						}
					}
				} else {
					// Ignore unmapped read and supplementary alignment
					if aread.Flags&sam.Unmapped != 0 || aread.Flags&sam.Supplementary != 0 {
						continue
					}
					pair.Reads = append(pair.Reads, aread)
				}
				sPair[iPair] = &pair
				if iPair == sPairLength-1 {
					select {
					case <-gctx.Done():
						return gctx.Err()
					case chAln <- sPair:
					}
					sPair = make([]*Pair, sPairLength)
					iPair = -1
				}
				iPair++
				nAlign++

				if verboseLevel > 0 {
					timeNow := time.Now()
					if timeNow.Sub(timeLog).Minutes() > 1. {
						fmt.Printf("%.1fmin - %s align. - %.2f Ma/hr\n", timeNow.Sub(timeStart).Minutes(), AddCommas(strconv.FormatUint(nAlign, 10)), (float64(nAlign)/timeNow.Sub(timeStart).Hours())/1000000.)
						timeLog = timeNow
					}
				}
			}
			// Send last packet
			if iPair > 0 {
				sPair = sPair[:iPair]
				select {
				case <-gctx.Done():
					return gctx.Err()
				case chAln <- sPair:
				}
			}
		}
		return nil
	})

	// Init cache pool
	pool := make(chan *Cache, nWorker2*2)
	for i := 0; i < cap(pool); i++ {
		c := NewCache(cacheLength, len(countMultis))
		pool <- c
	}

	// Spawn worker goroutine(s)
	g.Go(func() error {
		defer close(chFinal)
		// Start worker(s)
		wg, wgctx := errgroup.WithContext(gctx)
		for i := 0; i < nWorker2; i++ {
			wg.Go(func() error {
				var aread *sam.Record
				var apairKeep, coordProfileInside bool
				var pairCount float32
				var pairMulti int
				// Loop over data
				for sPair := range chAln {
					// Get cache
					c := <-pool
					for _, pair := range sPair {
						// Default to not keeping pair
						apairKeep = false

						// Alignment multiplicity
						if ignoreNHTag {
							pairMulti = 1
						} else {
							tag, found := pair.Reads[0].Tag([]byte{'N', 'H'})
							if found == false {
								return fmt.Errorf("Missing NH tag")
							} else {
								switch v := tag.Value().(type) {
								case uint8:
									pairMulti = int(v)
								case uint16:
									pairMulti = int(v)
								}
							}
						}
						pairCount = 1. / float32(pairMulti)

						// Input
						c.InputCount += 1. / float64(pairMulti)

						// Read length (both mates have to be desired length)
						if len(readLengths) > 0 {
							apairLengthOK := true
							for _, aread = range pair.Reads {
								areadLengthOK := false
								for _, l := range readLengths {
									if l == aread.Seq.Length {
										areadLengthOK = true
									}
								}
								if areadLengthOK == false {
									apairLengthOK = false
									break
								}
							}
							if apairLengthOK == false {
								continue
							}
						}

						// Read(s) filtering
						if inProperPair || minMappingQuality > 0 {
							filterOK := true
							for _, aread = range pair.Reads {
								// Is read in proper pair
								if inProperPair {
									if aread.Flags&sam.ProperPair == 0 {
										filterOK = false
										break
									}
								}
								// Minimum read mapping quality
								if minMappingQuality > 0 {
									if aread.MapQ < minMappingQuality {
										filterOK = false
										break
									}
								}
							}
							if filterOK == false {
								continue
							}
						}

						// Fragment length filtering
						if profileNoCoordMapping && (fragmentMinLength > 0 || fragmentMaxLength > 0) {
							fragmentLength := Abs(pair.Reads[0].TempLen)
							if fragmentMinLength > 0 && fragmentLength < fragmentMinLength {
								continue
							}
							if fragmentMaxLength > 0 && fragmentLength > fragmentMaxLength {
								continue
							}
						}

						// Read random selection
						if randProportion > 0. {
							if rand.Float32() > randProportion {
								continue
							}
						}

						// Get features overlap with reads
						featuresOverlap := feature.OverlapFeatureRead(pair.Reads, libraryR1Strand, trees)

						// Add reads to count and profile
						for featID, overlap := range featuresOverlap {
							if overlap.Length >= minOverlap {
								//if Debug {
								//	for i := 0; i < len(pair.Reads); i++ {
								//		alnRef, alnRead, alnSymbol := align.GetAln(pair.Reads[i])
								//		fmt.Println(string(alnRef))
								//		fmt.Println(string(alnSymbol))
								//		fmt.Println(string(alnRead), "\n")
								//	}
								//}
								// Feature
								feat := featureExts[featID]

								// Fragment length filtering
								if !profileNoCoordMapping && (fragmentMinLength > 0 || fragmentMaxLength > 0) {
									startProfile, endProfile := profile.FragmentCoords(pair.Reads, overlap, feat, profileNoCoordMapping)
									fragmentLength := endProfile - startProfile
									if fragmentMinLength > 0 && fragmentLength < fragmentMinLength {
										continue
									}
									if fragmentMaxLength > 0 && fragmentLength > fragmentMaxLength {
										continue
									}
								}

								// Increase cache size
								if len(c.Packets) <= c.LastPacket {
									c.Grow()
								}
								// Current feature
								c.Packets[c.LastPacket].ID = feat.ID

								// Profile
								if doProfile {
									coordProfileInside = false
									// Get read position within profile
									if pairMulti <= profileMulti {
										var err error
										switch profileType {
										case profile.ProfileTypeFirst:
											coordProfileInside, err = profile.ProfileFirst(pair.Reads, pair.OnlyRead1, paired, libraryR1Strand, overlap, feat, pairCount, c.Packets[c.LastPacket].ProfileChanges, profileNoCoordMapping, profileUntemplated, profileNoUntemplated)
										case profile.ProfileTypeLast:
											coordProfileInside = profile.ProfileLast(pair.Reads, pair.OnlyRead1, paired, libraryR1Strand, overlap, feat, pairCount, c.Packets[c.LastPacket].ProfileChanges, profileNoCoordMapping)
										case profile.ProfileTypeFirstLast:
											coordProfileInside = profile.ProfileFirstLast(pair.Reads, pair.OnlyRead1, paired, libraryR1Strand, overlap, feat, pairCount, c.Packets[c.LastPacket].ProfileChanges, profileNoCoordMapping)
										case profile.ProfileTypePosition:
											coordProfileInside = profile.ProfilePosition(pair.Reads, pair.OnlyRead1, paired, libraryR1Strand, overlap, feat, pairCount, c.Packets[c.LastPacket].ProfileChanges, profileNoCoordMapping, profilePositionFraction)
										case profile.ProfileTypeAll:
											coordProfileInside = profile.ProfileAll(pair.Reads, pair.OnlyRead1, paired, libraryR1Strand, overlap, feat, pairCount, c.Packets[c.LastPacket].ProfileChanges, profileNoCoordMapping)
										case profile.ProfileTypeSplice:
											coordProfileInside = profile.ProfileSplice(pair.Reads, pair.OnlyRead1, paired, libraryR1Strand, overlap, feat, pairCount, c.Packets[c.LastPacket].ProfileChanges, profileNoCoordMapping)
										case profile.ProfileTypeExtension:
											coordProfileInside = profile.ProfileExtension(pair.Reads, pair.OnlyRead1, paired, libraryR1Strand, overlap, feat, pairCount, c.Packets[c.LastPacket].ProfileChanges, profileNoCoordMapping, profileExtensionLength)
										}
										if err != nil {
											return err
										}
									}
									// Add read to profile
									if coordProfileInside {
										apairKeep = true
									}
								}

								// Count
								if countInProfile == false || coordProfileInside {
									for icm, cm := range countMultis {
										if pairMulti <= cm {
											apairKeep = true
											c.Packets[c.LastPacket].Counts[icm] += 1. / float64(pairMulti)
										}
									}
								}
								c.LastPacket++
							}
						}
						if apairKeep {
							iMulti := 0
							for icm, cm := range countMultis {
								if pairMulti <= cm {
									iMulti = icm
									break
								}
							}
							if countTotalRealRead {
								multiSets[iMulti].Add(aread.Name)
							} else {
								c.MultiCounts[iMulti] += 1. / float64(pairMulti)
							}
							if doOutSAM {
								for _, aread = range pair.Reads {
									samWriter.Write(aread)
								}
							}
						}
					}
					if c.LastPacket == 0 {
						pool <- c
					} else {
						select {
						case <-wgctx.Done():
							return wgctx.Err()
						case chFinal <- c:
						}
					}
				}
				return nil
			})
		}
		// Wait for the workers to finish
		err := wg.Wait()
		if err != nil {
			return err
		}
		return nil
	})

	// Combine data from worker into final count and profile
	nMulti := len(countMultis)
	for c := range chFinal {
		for i := 0; i < c.LastPacket; i++ {
			//DEBUG_PAIR fmt.Println("PACKET", i)
			// Count
			for j := 0; j < nMulti; j++ {
				featureExts[c.Packets[i].ID].Counts[1+(2*j)] += c.Packets[i].Counts[j]
				c.Packets[i].Counts[j] = 0
			}
			// Profile
			if doProfile {
				for j := 0; j <= c.Packets[i].ProfileChanges.ProfileLastIdx; j++ {
					//DEBUG_PAIR fmt.Println(c.Packets[i].ProfileChanges.ProfileIdxs[j], c.Packets[i].ID)
					featureExts[c.Packets[i].ID].Profile[c.Packets[i].ProfileChanges.ProfileIdxs[j]] += c.Packets[i].ProfileChanges.ProfileCounts[j]
				}
				c.Packets[i].ProfileChanges.ProfileLastIdx = -1
			}
		}
		// Total count
		for i := 0; i < nMulti; i++ {
			multisCounts[i] += c.MultiCounts[i]
			c.MultiCounts[i] = 0.
		}
		// Input count
		inputCount += c.InputCount
		c.InputCount = 0.
		// Reset
		c.LastPacket = 0
		pool <- c
	}

	err = g.Wait()
	if err != nil {
		return nAlign, err
	}

	// Normalization
	// Total length
	for i := 0; i < len(featureExts); i++ {
		countTotals[0] += featureExts[i].Counts[0]
	}
	// Total counts
	if !countTotalInput {
		var c int
		var p float64
		for icm := 0; icm < len(countMultis); icm++ {
			if countTotalRealRead {
				c += multiSets[icm].Size()
				countTotals[1+(2*icm)] = float64(c)
			} else {
				p += multisCounts[icm]
				countTotals[1+(2*icm)] = p
			}
		}
	}
	// Normalize counts to RPKM
	for icm := 0; icm < len(countMultis); icm++ {
		col := 1 + (2 * icm)
		if countTotals[col] > 0. {
			for _, feat := range featureExts {
				feat.Counts[col+1] = feat.Counts[col] * (1000. / feat.Counts[0]) * (1000000. / countTotals[col])
			}
			countTotals[col+1] = countTotals[col] * (1000. / countTotals[0]) * (1000000. / countTotals[col])
		}
	}
	// Normalize profiles to RPM
	if doProfile && profileNorm {
		normFactor := float32(1000000. / countTotals[profileMultiTotalCol])
		if verboseLevel > 0 {
			timeNow := time.Now()
			fmt.Printf("%.1fmin - Profile norm. factor: %f\n", timeNow.Sub(timeStart).Minutes(), normFactor)
		}
		for _, feat := range featureExts {
			for ip := 0; ip < len(feat.Profile); ip++ {
				feat.Profile[ip] *= normFactor
			}
		}
	}

	// Output: Count
	err = feature.WriteCounts(featureExts, countPath, countMultis, countTotals, appendOutput)
	if err != nil {
		return nAlign, err
	}
	// Output: Profile
	if doProfile {
		for ip := 0; ip < len(profileFormats); ip++ {
			if verboseLevel > 0 {
				timeNow := time.Now()
				fmt.Printf("%.1fmin - Writing %s output in %s\n", timeNow.Sub(timeStart).Minutes(), profileFormats[ip], profilePaths[ip])
			}
			feature.WriteProfiles(featureExts, featuresMapping, profilePaths[ip], profileFormats[ip], appendOutput)
		}
	}
	// Output: Report
	if pathReport != "" {
		err = WriteReport(pathReport, inputCount, countMultis, countTotalRealRead, multiSets, multisCounts)
		if err != nil {
			return nAlign, err
		}
	}

	return nAlign, nil
}
