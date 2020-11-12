// This file provides various utility functions needed by find-qubo.

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

// parseRow parses a row of Booleans, specified in a flexible manner, into a
// binary number.  It returns the number and the total number of bits on the
// line.  The function aborts on error.
func parseRow(s string) (uint, int) {
	// Discard comments ("#" to the end of the line).
	cIdx := strings.Index(s, "#")
	if cIdx != -1 {
		s = s[:cIdx]
	}

	// Ignore blank lines.
	fields := strings.Fields(s)
	n := len(fields)
	if n == 0 {
		return 0, 0
	}

	// Parse each field in turn.
	var v uint
	for _, f := range fields {
		v <<= 1
		switch strings.ToUpper(f) {
		case "0", "-1", "F", "FALSE":
		case "1", "T", "TRUE":
			v |= 1
		default:
			notify.Fatalf("Failed to parse %q as a Boolean value", f)
		}
	}
	return v, n
}

// ReadTruthTable reads a truth table from a file.
func ReadTruthTable(p *Parameters) TruthTable {
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
		v, nc := parseRow(ln)
		switch {
		case nc == 0:
			// Blank line: ignore.
			continue
		case prevNC == 0:
			// First row: Allocate the table, and store the first
			// value.
			prevNC = nc
			tt = make(TruthTable, 1<<nc)
			tt[v] = true
		case nc != prevNC:
			// Change in column count: Abort.
			notify.Fatalf("Column count changed from %d to %d in line %d", prevNC, nc, row)
		default:
			// Second or subsequent row: Store the value.
			tt[v] = true
		}
	}
	if err := scanner.Err(); err != nil {
		notify.Fatal(err)
	}
	return tt
}
