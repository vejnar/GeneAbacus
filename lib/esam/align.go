//
// Copyright (C) 2015-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package esam

import (
	"bytes"
	"strconv"
	"unicode"

	"github.com/biogo/hts/sam"
)

const (
	MDInsertion = iota
	MDMismatch
	MDSkip
)

type TagMDOp struct {
	Op     int
	Length int
	Seq    []byte
}

// ParseTagMD parses the MD attribute to blocks.
func ParseTagMD(rawTag string) (blocks []TagMDOp, err error) {
	var block []byte
	var l byte
	// Parsing tag
	i := 0
	for i < len(rawTag) {
		l = rawTag[i]
		if l == '^' {
			block = []byte("")
			i++ // Skipping "^"
			l = rawTag[i]
			for i < len(rawTag) {
				l = rawTag[i]
				if unicode.IsLetter(rune(l)) {
					block = append(block, l)
					i++
				} else {
					break
				}
			}
			blocks = append(blocks, TagMDOp{Op: MDInsertion, Length: len(block), Seq: block})
		} else if unicode.IsLetter(rune(l)) {
			blocks = append(blocks, TagMDOp{Op: MDMismatch, Length: 1, Seq: []byte{l}})
			i++
		} else {
			block = []byte("")
			for i < len(rawTag) {
				l = rawTag[i]
				if unicode.IsNumber(rune(l)) {
					block = append(block, l)
					i++
				} else {
					break
				}
			}
			step, err := strconv.Atoi(string(block))
			if err != nil {
				return blocks, nil
			}
			blocks = append(blocks, TagMDOp{Op: MDSkip, Length: step})
		}
	}
	return blocks, nil
}

// GetAln reconstitutes read alignment based on cigar string. The MD tag is used if present.
func GetAln(r *sam.Record) (ref, read, symbol []byte, err error) {
	// Parsing CIGAR string
	var iRead, length int
	var co sam.CigarOp
	var con sam.Consume
	seq := r.Seq.Expand()
	for i := 0; i < len(r.Cigar); i++ {
		co = r.Cigar[i]
		con = co.Type().Consumes()
		length = co.Len()
		if con.Query == 1 && con.Reference == 1 {
			ref = append(ref, seq[iRead:iRead+length]...)
			read = append(read, seq[iRead:iRead+length]...)
			if co.Type() == sam.CigarMatch {
				symbol = append(symbol, bytes.Repeat([]byte("|"), length)...)
			} else {
				symbol = append(symbol, bytes.Repeat([]byte("X"), length)...)
			}
			iRead += length
		} else if con.Query == 0 && con.Reference == 1 {
			ref = append(ref, bytes.Repeat([]byte("N"), length)...)
			read = append(read, bytes.Repeat([]byte("-"), length)...)
			symbol = append(symbol, bytes.Repeat([]byte("."), length)...)
		} else if con.Query == 1 && con.Reference == 0 {
			if co.Type() == sam.CigarInsertion {
				ref = append(ref, bytes.Repeat([]byte("-"), length)...)
				symbol = append(symbol, bytes.Repeat([]byte("."), length)...)
			} else {
				ref = append(ref, bytes.Repeat([]byte(" "), length)...)
				symbol = append(symbol, bytes.Repeat([]byte(" "), length)...)
			}
			read = append(read, seq[iRead:iRead+length]...)
			iRead += length
		}
	}
	// Parsing MD tag if present
	tag, found := r.Tag([]byte("MD"))
	if found {
		var blocks []TagMDOp
		blocks, err = ParseTagMD(tag.Value().(string))
		if err != nil {
			return
		}
		var iRef, iBlock int
		for _, b := range blocks {
			for iRef < len(ref) && (ref[iRef] == '-' || ref[iRef] == ' ') {
				iRef++
			}
			switch b.Op {
			case MDInsertion:
				for _, nt := range b.Seq {
					ref[iRef] = nt
					iRef++
				}
			case MDMismatch:
				ref[iRef] = b.Seq[0]
				symbol[iRef] = 'X'
				iRef++
			case MDSkip:
				iBlock = 0
				for iBlock < b.Length {
					if ref[iRef] != '-' && ref[iRef] != ' ' && symbol[iRef] != '.' {
						iBlock++
					}
					iRef++
				}
			}
		}
	}
	return
}
