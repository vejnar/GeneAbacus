//
// Copyright © 2015 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package feature

import (
	"bufio"
	"os"
	"strings"
)

func OpenMapping(mpath string) (map[string]string, error) {
	m := make(map[string]string)

	mfos, err := os.Open(mpath)
	if err != nil {
		return m, err
	}
	defer mfos.Close()

	tscanner := bufio.NewScanner(mfos)
	for tscanner.Scan() {
		fields := strings.Split(tscanner.Text(), "\t")
		m[fields[0]] = fields[1]
	}
	if err := tscanner.Err(); err != nil {
		return m, err
	}
	return m, nil
}

func MapName(name string, m map[string]string) string {
	if nn, ok := m[name]; ok {
		return nn
	}
	return name
}
