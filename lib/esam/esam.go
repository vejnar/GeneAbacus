//
// Copyright © 2015 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package esam

import (
	"github.com/biogo/hts/sam"
)

// PathSAM stores Path to SAM (Binary=false) or BAM (Binary=true) file.
type PathSAM struct {
	Path   string
	Binary bool
}

// Overlap returns the length of the overlap between the alignment of the SAM record and the interval specified with start and end.
func Overlap(r *sam.Record, start, end int) int {
	var overlap int
	pos := r.Pos
	for _, co := range r.Cigar {
		t := co.Type()
		con := t.Consumes()
		lr := co.Len() * con.Reference
		if con.Query == con.Reference {
			o := min(pos+lr, end) - max(pos, start)
			if o > 0 {
				overlap += o
			}
		}
		pos += lr
	}
	return overlap
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

// ShiftPos shifts the coordinate within the alignment. If strand is -1, it starts from the alignment end.
func ShiftPos(pos int, shift int, r *sam.Record, strand int8) (int, bool) {
	shifted := 0
	var con sam.Consume
	var lr int
	if strand == 1 {
		for i := 0; i < len(r.Cigar); i++ {
			co := r.Cigar[i]
			con = co.Type().Consumes()
			lr = co.Len() * con.Query
			if shifted+lr < shift {
				pos += co.Len() * con.Reference
				shifted += lr
			} else {
				if con.Reference != 0 {
					return pos + ((shift - shifted) * con.Reference), true
				} else {
					return 0, false
				}
			}
		}
	} else {
		for i := len(r.Cigar) - 1; i >= 0; i-- {
			co := r.Cigar[i]
			con = co.Type().Consumes()
			lr = co.Len() * con.Query
			if shifted+lr < shift {
				pos -= co.Len() * con.Reference
				shifted += lr
			} else {
				if con.Reference != 0 {
					return pos - ((shift - shifted) * con.Reference), true
				} else {
					return 0, false
				}
			}
		}
	}
	return 0, false
}
