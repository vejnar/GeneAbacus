//
// Copyright (C) 2015-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package profile

type ProfileChange struct {
	ProfileIdxs    []int
	ProfileCounts  []float32
	ProfileLastIdx int
}

func NewProfileChange(size int) *ProfileChange {
	c := ProfileChange{}
	c.ProfileLastIdx = -1
	c.ProfileIdxs = make([]int, size)
	c.ProfileCounts = make([]float32, size)
	return &c
}

func (c *ProfileChange) Write(i int, v float32) {
	c.ProfileLastIdx++
	if len(c.ProfileIdxs) <= c.ProfileLastIdx {
		c.Grow(2)
	}
	c.ProfileIdxs[c.ProfileLastIdx] = i
	c.ProfileCounts[c.ProfileLastIdx] = v
}

func (c *ProfileChange) Grow(factor int) {
	n := make([]int, len(c.ProfileIdxs)*factor)
	copy(n, c.ProfileIdxs)
	c.ProfileIdxs = n
	m := make([]float32, len(c.ProfileCounts)*factor)
	copy(m, c.ProfileCounts)
	c.ProfileCounts = m
}
