//
// Copyright (C) 2015-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//

package cmapper

type CoordMapper struct {
	CoordsGenome, CoordsTranscript [][]int
	Strand                         int8
	Length                         int
}

// Init.
func (cm *CoordMapper) Init() {
	// CoordsTranscript
	var tcoord int
	if cm.Strand == 1 {
		for i := 0; i < len(cm.CoordsGenome); i++ {
			exonLength := cm.CoordsGenome[i][1] - cm.CoordsGenome[i][0]
			cm.CoordsTranscript = append(cm.CoordsTranscript, []int{tcoord, tcoord + exonLength})
			tcoord += exonLength
		}
	} else if cm.Strand == -1 {
		for i := len(cm.CoordsGenome) - 1; i >= 0; i-- {
			exonLength := cm.CoordsGenome[i][1] - cm.CoordsGenome[i][0]
			cm.CoordsTranscript = append(cm.CoordsTranscript, []int{tcoord, tcoord + exonLength})
			tcoord += exonLength
		}
	}
	// Length
	cm.Length = cm.GetLength()
}

// GetLength returns mapper length.
func (cm *CoordMapper) GetLength() (length int) {
	for _, iv := range cm.CoordsGenome {
		length += iv[1] - iv[0]
	}
	return
}

// Genome2Transcript translates a coordinate from the genome to the transcript system.
func (cm *CoordMapper) Genome2Transcript(coord int) (tcoord int, within bool) {
	exonIth := -1
	for i := 0; i < len(cm.CoordsGenome); i++ {
		if coord >= cm.CoordsGenome[i][0] && coord < cm.CoordsGenome[i][1] {
			exonIth = i
			break
		}
	}
	if exonIth != -1 {
		if cm.Strand == 1 {
			tcoord = cm.CoordsTranscript[exonIth][0] + (coord - cm.CoordsGenome[exonIth][0])
		} else if cm.Strand == -1 {
			tcoord = cm.CoordsTranscript[len(cm.CoordsTranscript)-1-exonIth][1] - 1 - (coord - cm.CoordsGenome[exonIth][0])
		}
		within = true
	}
	return
}
