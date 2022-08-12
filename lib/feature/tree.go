//
// Copyright (C) 2015-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package feature

import (
	"github.com/biogo/hts/sam"
	"github.com/biogo/store/interval"

	"git.sr.ht/~vejnar/GeneAbacus/lib/esam"
)

type FeatureOverlap struct {
	Length int
	Read   []bool
}

// BuildFeatTrees builds a tree of features: each interval (i.e. exon) of each feature is added to the tree.
func BuildFeatTrees(features []Feature) (trees map[string]map[int8]*interval.IntTree, err error) {
	trees = make(map[string]map[int8]*interval.IntTree)
	icoord := 0
	for _, feat := range features {
		for _, coord := range feat.Coords {
			// New tree for unseen chromosome
			if _, ok := trees[feat.Chrom]; !ok {
				trees[feat.Chrom] = make(map[int8]*interval.IntTree)
				trees[feat.Chrom][1] = &interval.IntTree{}
				trees[feat.Chrom][-1] = &interval.IntTree{}
			}
			// Creating new interval
			iv := IntInterval{Start: coord[0], End: coord[1], UID: uintptr(icoord), Feature: Feature{ID: feat.ID, Name: feat.Name, Strand: feat.Strand}}
			// Inserting interval
			err = trees[feat.Chrom][feat.Strand].Insert(iv, false)
			if err != nil {
				return
			}
			icoord++
		}
	}
	for k, _ := range trees {
		trees[k][1].AdjustRanges()
		trees[k][-1].AdjustRanges()
	}
	return
}

func OverlapFeatureRead(areads []*sam.Record, libraryR1Strand int8, trees map[string]map[int8]*interval.IntTree) map[uint32]FeatureOverlap {
	// Read strand corrected for library strand
	var areadStrands []int8
	apairR1Strand := areads[0].Strand()
	if libraryR1Strand == 0 {
		areadStrands = []int8{-1, 1}
	} else if libraryR1Strand == 1 {
		areadStrands = []int8{apairR1Strand}
	} else if libraryR1Strand == -1 {
		areadStrands = []int8{apairR1Strand * -1}
	}
	// Aligned-read intervals
	var areadIntervals []IntInterval
	for _, aread := range areads {
		areadIntervals = append(areadIntervals, IntInterval{Start: aread.Start(), End: aread.End()})
	}
	// Overlapping features
	featuresOverlap := make(map[uint32]FeatureOverlap)
	for i := 0; i < len(areads); i++ {
		for _, rstrand := range areadStrands {
			if tree, ok := trees[areads[i].Ref.Name()]; ok {
				for _, iv := range tree[rstrand].Get(areadIntervals[i]) {
					if fo, ok := featuresOverlap[iv.(IntInterval).Feature.ID]; ok {
						fo.Length += esam.Overlap(areads[i], iv.Range().Start, iv.Range().End)
						fo.Read[i] = true
					} else {
						featuresOverlap[iv.(IntInterval).Feature.ID] = FeatureOverlap{Length: esam.Overlap(areads[i], iv.Range().Start, iv.Range().End), Read: make([]bool, len(areads))}
						featuresOverlap[iv.(IntInterval).Feature.ID].Read[i] = true
					}
					//if Debug {
					//	fmt.Println("Overlap", "read"+strconv.Itoa(i+1), iv.(IntInterval).Feature.Name, iv.(IntInterval).UID, iv.Range().Start, iv.Range().End, Overlap(areads[i], iv.Range().Start, iv.Range().End))
					//}
				}
			}
		}
	}
	return featuresOverlap
}
