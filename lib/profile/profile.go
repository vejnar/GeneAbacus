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

	"git.sr.ht/~vejnar/GeneAbacus/lib/esam"
)

const (
	ProfileTypeNone = iota
	ProfileTypeFirst
	ProfileTypeLast
	ProfileTypeFirstLast
	ProfileTypePosition
	ProfileTypeAll
	ProfileTypeSplice
	ProfileTypeExtension
)

// TrimUntemplated returns the number of untemplated nucleotide (max of maxShift).
func TrimUntemplated(r *sam.Record, maxShift int, strand int8) (int, error) {
	var iSymbol, iMatch, iMismatch, lenTU int
	iMismatch = -1
	_, _, symbol, err := esam.GetAln(r)
	if err != nil {
		return lenTU, err
	}
	if strand == 1 {
		for iMatch < maxShift && iSymbol < len(symbol) {
			if symbol[iSymbol] != '.' {
				if symbol[iSymbol] != '|' {
					iMismatch = iMatch
				}
				iMatch++
			}
			iSymbol++
		}
		lenTU = iMismatch + 1
	} else {
		iSymbol = len(symbol) - 1
		for iMatch < maxShift && iSymbol > 0 {
			if symbol[iSymbol] != '.' {
				if symbol[iSymbol] != '|' {
					iMismatch = iMatch
				}
				iMatch++
			}
			iSymbol--
		}
		lenTU = iMismatch + 1
	}
	return lenTU, err
}
