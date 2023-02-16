//
// Copyright © 2015 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package profile

import (
	"math"

	"github.com/biogo/hts/sam"

	"git.sr.ht/~vejnar/GeneAbacus/lib/feature"
)

func FragmentCoords(areads []*sam.Record, overlap feature.FeatureOverlap, feat *feature.FeatureExt, profileNoCoordMapping bool) (int, int) {
	// Search start and end of alignment on feature
	var startProfile, endProfile, coordProfile int
	startProfileGenome := math.MaxInt32
	endProfileGenome := math.MinInt32
	coordProfileInside := true
	// Find first mapped position
	for i := 0; i < len(areads); i++ {
		if overlap.Read[i] {
			for iar := areads[i].Start(); iar < areads[i].End(); iar++ {
				if profileNoCoordMapping {
					coordProfile = iar
				} else {
					coordProfile, coordProfileInside = feat.CoordMapper.Genome2Transcript(iar)
				}
				if coordProfileInside && iar < startProfileGenome {
					startProfileGenome = iar
					startProfile = coordProfile
					break
				}
			}
		}
	}
	// Find last mapped position
	for i := 0; i < len(areads); i++ {
		if overlap.Read[i] {
			for iar := areads[i].End() - 1; iar >= areads[i].Start(); iar-- {
				if profileNoCoordMapping {
					coordProfile = iar
				} else {
					coordProfile, coordProfileInside = feat.CoordMapper.Genome2Transcript(iar)
				}
				if coordProfileInside && iar > endProfileGenome {
					endProfileGenome = iar
					endProfile = coordProfile
					break
				}
			}
		}
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
	return startProfile, endProfile
}

func ProfileAll(areads []*sam.Record, onlyRead1 bool, paired bool, libraryR1Strand int8, overlap feature.FeatureOverlap, feat *feature.FeatureExt, pairCount float32, profileChanges *ProfileChange, profileNoCoordMapping bool) bool {
	// Get start and end coordinates of fragment
	startProfile, endProfile := FragmentCoords(areads, overlap, feat, profileNoCoordMapping)
	// Add count
	for ip := startProfile; ip < endProfile; ip++ {
		profileChanges.Write(ip, pairCount)
	}
	return true
}
