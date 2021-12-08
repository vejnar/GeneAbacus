//
// Copyright (C) 2015-2021 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package feature

import (
	"fmt"

	"github.com/biogo/store/interval"
)

// Integer-specific intervals

type IntInterval struct {
	Start, End int
	UID        uintptr
	Feature    Feature
}

func (i IntInterval) OverlapLength(b interval.IntRange) int {
	return min(i.End, b.End) - max(i.Start, b.Start)
}

func (i IntInterval) Overlap(b interval.IntRange) bool {
	// Half-open interval indexing.
	return i.End > b.Start && i.Start < b.End
}

func (i IntInterval) ID() uintptr {
	return i.UID
}

func (i IntInterval) Range() interval.IntRange {
	return interval.IntRange{i.Start, i.End}
}

func (i IntInterval) String() string {
	return fmt.Sprintf("[%d,%d)#%d-%s", i.Start, i.End, i.UID, i.Feature.Name)
}

func min(a, b int) int {
	if a > b {
		return b
	}
	return a
}

func max(a, b int) int {
	if a < b {
		return b
	}
	return a
}
