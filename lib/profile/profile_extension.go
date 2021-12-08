//
// Copyright (C) 2015-2021 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package profile

import (
	"github.com/biogo/hts/sam"

	"github.com/vejnar/geneabacus/lib/feature"
)

// ProfileExtension is designed for unstranded single-end sequencing, so read2 and libraryR1Strand are ignored
func ProfileExtension(areads []*sam.Record, onlyRead1 bool, paired bool, libraryR1Strand int8, overlap feature.FeatureOverlap, feat *feature.FeatureExt, pairCount float32, profileChanges *ProfileChange, profileNoCoordMapping bool, profileExtensionLength int) bool {
	var startProfile, endProfile int
	var coordProfileInside bool
	// Check that overlap is for this read
	if !overlap.Read[0] {
		return false
	}
	// Get genomic position
	if areads[0].Strand() == 1 {
		startProfile = areads[0].Start()
		// Coordinate must not go beyong feature end
		endProfile = min(feat.Coords[len(feat.Coords)-1][1], areads[0].Start()+profileExtensionLength)
	} else {
		startProfile = max(0, areads[0].End()-profileExtensionLength)
		endProfile = areads[0].End()
	}
	if !profileNoCoordMapping {
		var startProfileInside, endProfileInside bool
		// Transpose from genome to transcript coordinate
		startProfile, startProfileInside = feat.CoordMapper.Genome2Transcript(startProfile)
		endProfile, endProfileInside = feat.CoordMapper.Genome2Transcript(endProfile - 1)
		if startProfileInside || endProfileInside {
			coordProfileInside = true
		}
		// Check orientation
		if startProfile < endProfile {
			// Half-open end
			endProfile++
		} else {
			// Minus strand and half-open end
			s := startProfile
			startProfile = endProfile
			endProfile = s + 1
		}
	}
	// Add count
	for ip := startProfile; ip < endProfile; ip++ {
		profileChanges.Write(ip, pairCount)
	}
	return coordProfileInside
}

func min(a, b int) int {
	if a > b {
		return b
	}
	return a
}
