// This file provides functions for manipulating truth tables.

package main

import (
	"bufio"
	"io"
	"os"
	"strings"
)

// A TruthTable represents a complete truth table with each row marked as
// either valid or invalid.
type TruthTable []bool

// parseRow parses a row of Booleans, specified in a flexible manner, into one
// or more binary numbers.  It returns the numbers and the total number of bits
// on the line.  The function aborts on error.
func parseRow(s string) ([]uint, int) {
	// Discard comments ("#" to the end of the line).
	cIdx := strings.Index(s, "#")
	if cIdx != -1 {
		s = s[:cIdx]
	}

	// Ignore blank lines.
	fields := strings.Fields(s)
	n := len(fields)
	if n == 0 {
		return nil, 0
	}

	// Parse each field in turn.
	vs := make([]uint, 1, 128)
	for _, f := range fields {
		for i := range vs {
			vs[i] <<= 1
		}
		switch strings.ToUpper(f) {
		case "0", "-1", "F", "FALSE":
			// False
		case "1", "T", "TRUE":
			// True
			for i := range vs {
				vs[i] |= 1
			}
		case "*", "-", "?":
			// Don't care
			for i := range vs {
				vs = append(vs, vs[i]|1)
			}
		default:
			notify.Fatalf("Failed to parse %q as a Boolean value", f)
		}
	}
	return vs, n
}

// ReadTruthTable reads a truth table from a file.  It returns both the truth
// table and a column count.  This function aborts on error.
func ReadTruthTable(p *Parameters) (TruthTable, int) {
	// Open the input file.
	var r io.Reader
	if p.TTName == "" {
		// Read from standard input.
		r = os.Stdin
	} else {
		// Read from the named file.
		f, err := os.Open(p.TTName)
		if err != nil {
			notify.Fatal(err)
		}
		defer f.Close()
		r = f
	}

	// Read and parse each line in turn.
	var tt TruthTable
	var prevNC int
	row := 0
	scanner := bufio.NewScanner(r)
	for scanner.Scan() {
		row++
		ln := scanner.Text()
		vs, nc := parseRow(ln)
		switch {
		case nc == 0:
			// Blank line: ignore.
			continue
		case prevNC == 0:
			// First row: Allocate the table, and store the first
			// value.
			prevNC = nc
			tt = make(TruthTable, 1<<nc)
			for _, v := range vs {
				tt[v] = true
			}
		case nc != prevNC:
			// Change in column count: Abort.
			notify.Fatalf("Column count changed from %d to %d in line %d", prevNC, nc, row)
		default:
			// Second or subsequent row: Store the value.
			for _, v := range vs {
				tt[v] = true
			}
		}
	}
	if err := scanner.Err(); err != nil {
		notify.Fatal(err)
	}
	if prevNC == 0 {
		notify.Fatal("Truth table is empty")
	}

	// Append ancilla columns to the table we just created.
	return tt.AppendAncillae(prevNC, p.NAnc), prevNC + p.NAnc
}

// AppendAncillae returns a new truth table with ancillary columns appended to
// the right.
func (tt TruthTable) AppendAncillae(nc, na int) TruthTable {
	// Return the original truth table if we have no ancillae to add.
	if na == 0 {
		return tt
	}

	// Iterate over each row of the original truth table.
	tta := make(TruthTable, 0, 1<<(nc+na))
	for _, b := range tt {
		// Replicate each row 2^na times.
		for i := 0; i < 1<<na; i++ {
			tta = append(tta, b)
		}
	}
	return tta
}
