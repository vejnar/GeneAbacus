//
// Copyright (C) 2015-2021 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package main

import (
	"encoding/json"
	"fmt"
	"os"

	"gopkg.in/fatih/set.v0"
)

func WriteReport(pathReport string, countMultis []int, countTotalRealRead bool, multiSets []set.Interface, multisCounts []float64) (err error) {
	countReport := make(map[string]uint32)
	for i := 0; i < len(countMultis); i++ {
		if countMultis[i] == 1 {
			if countTotalRealRead {
				countReport["align_unique"] = uint32(multiSets[i].Size())
			} else {
				countReport["align_unique"] = uint32(multisCounts[i])
			}
		} else {
			if countTotalRealRead {
				countReport["align_multi"] += uint32(multiSets[i].Size())
			} else {
				countReport["align_multi"] += uint32(multisCounts[i])
			}
		}
	}
	report, _ := json.MarshalIndent(countReport, "", "  ")
	if pathReport != "-" {
		if f, err := os.Create(pathReport); err != nil {
			return err
		} else {
			f.Write(report)
			f.Close()
		}
	} else {
		fmt.Println(string(report))
	}
	return nil
}
