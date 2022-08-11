//
// Copyright (C) 2015-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package profile

import (
	"math"

	"github.com/biogo/hts/sam"

	"github.com/vejnar/geneabacus/lib/feature"
)

func ProfilePosition(areads []*sam.Record, onlyRead1 bool, paired bool, libraryR1Strand int8, overlap feature.FeatureOverlap, feat *feature.FeatureExt, pairCount float32, profileChanges *ProfileChange, profileNoCoordMapping bool, profilePositionFraction float64) bool {
	// Get start and end coordinates of fragment
	startProfile, endProfile := FragmentCoords(areads, overlap, feat, profileNoCoordMapping)
	// Add count
	profileChanges.Write(startProfile+int(math.Round(float64(endProfile-1-startProfile)*profilePositionFraction)), pairCount)
	return true
}
