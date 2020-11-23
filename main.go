/*
Find parameters for a QUBO given a truth table.
*/

package main

import (
	"fmt"
	"log"
	"os"
	"sort"
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
	notify = log.New(os.Stderr, os.Args[0]+": ", 0)
	status = log.New(os.Stderr, "INFO: ", 0)
	var p Parameters
	ParseCommandLine(&p)
	p.SeparatedGen = -1 // We haven't yet separated valid from invalid.

	// Read the input file.
	tt, nc := ReadTruthTable(&p)
	p.TT = tt
	p.NCols = nc

	// Precompute a matrix with all possible 0/1 columns.
	p.AllCols = AllPossibleColumns(p.NCols)

	// Precompute the maximum gap, rounded up to a power of 10.
	p.MaxGap = findMaxGap(&p)

	// Try to find coefficients that represent the truth table.
	qubo, bad, nGen := OptimizeCoeffs(&p)

	// Output what we found.
	fmt.Printf("Total number of generations = %d\n", nGen)
	fmt.Printf("Final badness = %v\n", bad)
	fmt.Printf("Final coefficients = %v\n", qubo.Coeffs)
	qubo.Rescale()
	fmt.Printf("Rescaled coefficients = %v\n", qubo.Coeffs)
	qubo.Evaluate() // Recompute the gap.
	fmt.Printf("Final valid/invalid gap = %v\n", qubo.Gap)
	fmt.Printf("Matrix form = %v\n", qubo.AsOctaveMatrix())
	vals := qubo.EvaluateAllInputs()
	isValid := qubo.SelectValidRows(vals)
	outputEvaluation(&p, isValid, vals)
}
