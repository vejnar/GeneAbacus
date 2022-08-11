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

	"github.com/vejnar/geneabacus/lib/feature"
)

func ProfileSplice(areads []*sam.Record, onlyRead1 bool, paired bool, libraryR1Strand int8, overlap feature.FeatureOverlap, feat *feature.FeatureExt, pairCount float32, profileChanges *ProfileChange, profileNoCoordMapping bool) bool {
	var length, coordProfile int
	var co sam.CigarOp
	var con sam.Consume
	var coordProfileInside, inside, newCoord bool
	profileLastIdxStart := max(0, profileChanges.ProfileLastIdx)
	for ir := 0; ir < len(areads); ir++ {
		if overlap.Read[ir] {
			iRef := areads[ir].Start()
			for ic := 0; ic < len(areads[ir].Cigar); ic++ {
				co = areads[ir].Cigar[ic]
				con = co.Type().Consumes()
				length = co.Len()
				if con.Query == 1 && con.Reference == 1 {
					if profileNoCoordMapping {
						coordProfileInside = true
					}
					// Add read part
					for ip := iRef; ip < iRef+length; ip++ {
						newCoord = true
						if profileNoCoordMapping {
							coordProfile = ip
						} else {
							coordProfile, inside = feat.CoordMapper.Genome2Transcript(ip)
							if inside && coordProfileInside == false {
								coordProfileInside = true
							}
						}
						for id := profileLastIdxStart; id <= profileLastIdxStart+profileChanges.ProfileLastIdx; id++ {
							if profileChanges.ProfileIdxs[id] == coordProfile {
								newCoord = false
							}
						}
						if newCoord {
							profileChanges.Write(coordProfile, pairCount)
						}
					}
				}
				if con.Reference == 1 {
					iRef += length
				}
			}
		}
	}
	return coordProfileInside
}

func max(a, b int) int {
	if a < b {
		return b
	}
	return a
}
