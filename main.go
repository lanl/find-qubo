/*
Find parameters for a QUBO given a truth table.
*/

package main

import (
	"fmt"
	"log"
	"os"
	"sort"
	"time"
)

// notify is used to output error messages.
var notify *log.Logger

// status is used to output status messages.
var status *log.Logger

// outputEvaluation pretty-prints the evaluation of all inputs.
func outputEvaluation(p *Parameters, isValid []bool, eval []float64) {
	// Map each value to its rank.
	sorted := make([]float64, len(eval))
	copy(sorted, eval)
	sort.Float64s(sorted)
	rank := make(map[float64]int)
	for i, v := range sorted {
		if _, ok := rank[v]; ok {
			continue // Ignore duplicate values
		}
		rank[v] = i + 1
	}

	// Tally the number of valid rows.
	nValid := 0
	for _, v := range isValid {
		if v {
			nValid++
		}
	}

	// Output each input string, output value, and rank.
	status.Print("Complete evaluation:")
	digits := len(fmt.Sprintf("%d", len(eval)+1))
	for i, v := range eval {
		// Set validMark to "*" for valid rows, " " for invalid rows.
		validMark := ' '
		if isValid[i] {
			validMark = '*'
		}

		// Set badRank to "X" for misordered ranks, " " for correct
		// ranks.
		badRank := ' '
		switch {
		case isValid[i] && rank[v] > nValid:
			badRank = 'X'
		case !isValid[i] && rank[v] <= nValid:
			badRank = 'X'
		}

		// Output the current row of the truth table.
		fmt.Printf("    %0*b %c  %18.15f  %*d %c\n", p.NCols, i, validMark, v, digits, rank[v], badRank)
	}
}

func main() {
	// Initialize program parameters.
	startTime := time.Now()
	notify = log.New(os.Stderr, os.Args[0]+": ", 0)
	status = log.New(os.Stderr, "INFO: ", 0)
	var p Parameters
	ParseCommandLine(&p)

	// Read the input file.
	p.TT, p.NCols = ReadTruthTable(&p)

	// Precompute a matrix with all possible 0/1 columns.
	p.AllCols = AllPossibleColumns(p.NCols)

	// Find the best QUBO we can.
	q, vals, isValid := OptimizeCoeffs(&p)
	if q == nil {
		// TODO: Increase the number of ancillae and try again.
		notify.Fatal("No solution was found.")
	}

	// Output what we found.
	status.Printf("Total program run time: %v", time.Since(startTime))
	fmt.Printf("Final coefficients = %v\n", q.Coeffs)
	fmt.Printf("Matrix form = %v\n", q.AsOctaveMatrix())
	outputEvaluation(&p, isValid, vals)
}
