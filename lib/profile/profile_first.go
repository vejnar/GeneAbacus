//
// Copyright (C) 2015-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package profile

import (
	"github.com/biogo/hts/sam"

	"github.com/vejnar/geneabacus/lib/esam"
	"github.com/vejnar/geneabacus/lib/feature"
)

func ProfileFirst(areads []*sam.Record, onlyRead1 bool, paired bool, libraryR1Strand int8, overlap feature.FeatureOverlap, feat *feature.FeatureExt, pairCount float32, profileChanges *ProfileChange, profileNoCoordMapping bool, profileUntemplated int, profileNoUntemplated bool) (bool, error) {
	var coordProfile, iRead int
	var coordProfileInside bool
	// Determine which read is first in case of paired-end sequencing
	if paired {
		iRead = -1
		if libraryR1Strand == 1 && (onlyRead1 || len(areads) == 2) {
			iRead = 0
		} else if libraryR1Strand == -1 {
			if len(areads) == 2 {
				iRead = 1
			} else if len(areads) == 1 && !onlyRead1 {
				iRead = 0
			}
		}
	}
	// Compute where to add the read
	if iRead != -1 {
		// Check that overlap is for this read
		if !overlap.Read[iRead] {
			return coordProfileInside, nil
		}
		// Get genomic position
		if feat.Strand == 1 {
			coordProfile = areads[iRead].Start()
		} else {
			coordProfile = areads[iRead].End() - 1
		}
		// Untemplated nucleotide
		if profileUntemplated > 0 {
			lenTU, err := TrimUntemplated(areads[iRead], profileUntemplated, feat.Strand)
			if err != nil {
				return coordProfileInside, err
			}
			//if Debug {
			//	fmt.Println(areads[iRead].Name)
			//	for i := 0; i < len(areads); i++ {
			//		alnRef, alnRead, alnSymbol := esam.GetAln(areads[i])
			//		fmt.Println(string(alnRef))
			//		fmt.Println(string(alnSymbol))
			//		fmt.Println(string(alnRead))
			//	}
			//	fmt.Println("Trimming untemplated:", areads[iRead].Name, lenTU, coordProfile, feat.Strand)
			//}
			if lenTU > 0 {
				if profileNoUntemplated {
					coordProfileInside = false
				} else {
					coordProfile, coordProfileInside = esam.ShiftPos(coordProfile, lenTU, areads[iRead], feat.Strand)
				}
			} else {
				coordProfileInside = true
			}
			//if Debug {
			//	fmt.Println("Trimming untemplated, New coordinate:", coordProfile, "\n")
			//}
		} else {
			coordProfileInside = true
		}
		if coordProfileInside && !profileNoCoordMapping {
			// Transpose from genome to transcript coordinate
			coordProfile, coordProfileInside = feat.CoordMapper.Genome2Transcript(coordProfile)
		}
		// Add count
		if coordProfileInside {
			profileChanges.Write(coordProfile, pairCount)
		}
	}
	return coordProfileInside, nil
}
