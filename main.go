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

	// Set the remaining parameters.
	PrepareGAParameters(&p)

	// Try to find coefficients that represent the truth table.
	var qubo *QUBO  // Best QUBO found
	var totGen uint // Number of generations evolved across all GAs run
	for {
		q, bad, nGen := OptimizeCoeffs(&p)
		totGen += nGen
		if q.Gap > 0.0 {
			// We found a valid solution!
			status.Print("We found a valid, but likely suboptimal, solution!")
			status.Printf("Final badness = %v", bad)
			status.Printf("Total generations across all executed GAs = %d", totGen)
			qubo = q
			break
		}

		// We failed to separate valid from invalid rows.  See if
		// adding an ancillary variable helps.
		varStr := "variables"
		if p.NAnc == 1 {
			varStr = "variable"
		}
		status.Printf("A solution with %d ancillary %s seems unlikely.", p.NAnc, varStr)
		status.Printf("Increasing the number of ancillae from %d to %d and restarting the genetic algorithm.", p.NAnc, p.NAnc+1)
		p.TT = p.TT.AppendAncillae(p.NCols, 1)
		p.NCols++
		p.NAnc++
		PrepareGAParameters(&p)
	}

	// Report GA statistics.
	if len(qubo.History) > 0 {
		status.Print("GA mutation statistics:")
		tot := 0
		for _, v := range qubo.History {
			tot += v
		}
		for k, v := range qubo.History {
			status.Printf("    %-25s %10d (%5.1f%%)", k, v, 100.0*float64(v)/float64(tot))
		}
	}

	// Now that we've separated valid from invalid rows, use a
	// linear-programming solver to ensure valid-row equality and to
	// maximize the invalid-valid gap.
	status.Print("Using linear programming to improve the solution.")
	vals := qubo.EvaluateAllInputs()
	isValid := qubo.SelectValidRows(vals)
	ok := qubo.LPReoptimize(isValid)
	if ok {
		qubo.Evaluate()                 // Recompute the gap.
		vals = qubo.EvaluateAllInputs() // Recompute each row's value.
	} else {
		status.Print("The LP solver failed to improve the solution.")
		qubo.Rescale() // Not needed if LP succeeded.
	}

	// Output what we found.
	status.Printf("Total program run time: %v", time.Since(startTime))
	fmt.Printf("Final coefficients = %v\n", qubo.Coeffs)
	fmt.Printf("Final valid/invalid gap = %v\n", qubo.Gap)
	fmt.Printf("Matrix form = %v\n", qubo.AsOctaveMatrix())
	outputEvaluation(&p, isValid, vals)
}
