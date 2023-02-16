//
// Copyright © 2015 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package profile

import (
	"github.com/biogo/hts/sam"

	"git.sr.ht/~vejnar/GeneAbacus/lib/feature"
)

func ProfileLast(areads []*sam.Record, onlyRead1 bool, paired bool, libraryR1Strand int8, overlap feature.FeatureOverlap, feat *feature.FeatureExt, pairCount float32, profileChanges *ProfileChange, profileNoCoordMapping bool) bool {
	var coordProfile, iRead int
	var coordProfileInside bool
	// Determine which read is last in case of paired-end sequencing
	if paired {
		iRead = -1
		if libraryR1Strand == 1 {
			if len(areads) == 1 && !onlyRead1 {
				iRead = 0
			} else if len(areads) == 2 {
				iRead = 1
			}
		} else if libraryR1Strand == -1 && (onlyRead1 || len(areads) == 2) {
			iRead = 0
		}
	}
	// Compute where to add the read
	if iRead != -1 {
		// Check that overlap is for this read
		if !overlap.Read[iRead] {
			return coordProfileInside
		}
		// Get genomic position
		if feat.Strand == 1 {
			coordProfile = areads[iRead].End() - 1
		} else {
			coordProfile = areads[iRead].Start()
		}
		// Transpose from genome to transcript coordinate
		if profileNoCoordMapping {
			coordProfileInside = true
		} else {
			coordProfile, coordProfileInside = feat.CoordMapper.Genome2Transcript(coordProfile)
		}
		// Add count
		if coordProfileInside {
			profileChanges.Write(coordProfile, pairCount)
		}
	}
	return coordProfileInside
}
