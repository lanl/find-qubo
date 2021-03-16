// This file provides functions for manipulating truth tables.

package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"
)

// A TruthTable represents a complete truth table (2^n rows for n columns) with
// each row marked as either valid or invalid.
type TruthTable struct {
	TT    []bool // The truth table proper, with each row representing "valid" or "invalid"
	NCols int    // Number of columns
	NRows int    // Number of rows
}

// NewTruthTable returns an empty truth table (i.e., all rows are marked
// invalid).
func NewTruthTable(nc int) TruthTable {
	var tt TruthTable
	tt.NCols = nc
	tt.NRows = 1 << nc
	tt.TT = make([]bool, tt.NRows)
	return tt
}

// Copy deep-copies a given truth table.
func (tt TruthTable) Copy() TruthTable {
	ttc := tt // Start by copying all fields.
	ttc.TT = make([]bool, len(tt.TT))
	copy(ttc.TT, tt.TT)
	return ttc
}

// Clear marks all truth-table rows as invalid.
func (tt TruthTable) Clear() {
	for i := range tt.TT {
		tt.TT[i] = false
	}
}

// parseRow parses a single row of Booleans, specified in a flexible manner,
// into a binary numbers.  It returns both the number and the total number of
// bits on the line.
func parseRow(s string) (uint, int, error) {
	// Discard comments ("#" to the end of the line).
	cIdx := strings.Index(s, "#")
	if cIdx != -1 {
		s = s[:cIdx]
	}

	// Ignore blank lines.
	fields := strings.Fields(s)
	n := len(fields)
	if n == 0 {
		return 0, 0, nil
	}

	// Parse each field in turn.
	var v uint
	for _, f := range fields {
		v <<= 1
		switch strings.ToUpper(f) {
		case "0", "-1", "F", "FALSE":
			// False
		case "1", "T", "TRUE":
			// True
			v |= 1
		default:
			return 0, 0, fmt.Errorf("failed to parse %q as a Boolean value", f)
		}
	}
	return v, n, nil
}

// ReadTruthTable reads a truth table from a file.
func ReadTruthTable(fn string) (TruthTable, error) {
	// Open the input file.
	var r io.Reader
	if fn == "" {
		// Read from standard input.
		r = os.Stdin
	} else {
		// Read from the named file.
		f, err := os.Open(fn)
		if err != nil {
			return TruthTable{}, err
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
		v, nc, err := parseRow(ln)
		if err != nil {
			return TruthTable{}, err
		}
		switch {
		case nc == 0:
			// Blank line: ignore.
			continue
		case prevNC == 0:
			// First row: Allocate the table, and store the first
			// value.
			prevNC = nc
			tt = NewTruthTable(nc)
			tt.TT[v] = true
		case nc != prevNC:
			// Change in column count: Abort.
			return TruthTable{}, fmt.Errorf("column count changed from %d to %d in line %d", prevNC, nc, row)
		default:
			// Second or subsequent row: Store the value.
			tt.TT[v] = true
		}
	}
	if err := scanner.Err(); err != nil {
		return TruthTable{}, err
	}
	if prevNC == 0 {
		return TruthTable{}, fmt.Errorf("truth table is empty")
	}
	return tt, nil
}

// ValidRows returns a list of indices of all valid rows.
func (tt *TruthTable) ValidRows() []int {
	rows := make([]int, 0, tt.NRows)
	for i, v := range tt.TT {
		if v {
			rows = append(rows, i)
		}
	}
	return rows
}
