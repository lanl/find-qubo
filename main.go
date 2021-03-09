/*
Find parameters for a QUBO given a truth table.
*/

package main

import (
	"fmt"
	"log"
	"os"
	"runtime/pprof"
	"sort"
	"strings"
)

// notify is used to output error messages.
var notify *log.Logger

// info is used to output status messages.
var info *log.Logger

// outputEvaluation pretty-prints the evaluation of all inputs.
func outputEvaluation(p *Parameters, tt TruthTable, eval []float64) {
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
	for _, v := range tt.TT {
		if v {
			nValid++
		}
	}

	// Output each input string, output value, and rank.
	fmt.Println("Complete evaluation:")
	digits := len(fmt.Sprintf("%d", len(eval)+1))
	for i, v := range eval {
		// Set validMark to "*" for valid rows, " " for invalid rows.
		validMark := ' '
		if tt.TT[i] {
			validMark = '*'
		}

		// Set badRank to "X" for misordered ranks, " " for correct
		// ranks.
		badRank := ' '
		switch {
		case tt.TT[i] && rank[v] > nValid:
			badRank = 'X'
		case !tt.TT[i] && rank[v] <= nValid:
			badRank = 'X'
		}

		// Output the current row of the truth table.
		fmt.Printf("    %0*b %c  %18.15f  %*d %c\n", tt.NCols, i, validMark, v, digits, rank[v], badRank)
	}
}

// Return a list of coefficients as GNU Octave or MATLAB matrix input.
func coeffsToOctave(n int, coeffs []float64) string {
	oct := make([]string, n)
	i := n
	for r := 0; r < n; r++ {
		row := make([]string, n)
		for c := 0; c < n; c++ {
			switch {
			case c < r:
				// Not in upper triangle
				row[c] = "0"
			case c == r:
				// Linear term
				row[c] = fmt.Sprint(coeffs[c])
			default:
				// Quadratic term
				row[c] = fmt.Sprint(coeffs[i])
				i++
			}
		}
		oct[r] = strings.Join(row, " ")
	}
	return "[" + strings.Join(oct, " ; ") + "]"
}

func main() {
	// Initialize program parameters.
	notify = log.New(os.Stderr, os.Args[0]+": ", 0)
	info = log.New(os.Stderr, "INFO: ", 0)
	var p Parameters
	ParseCommandLine(&p)

	// Read the input file.
	tt, err := ReadTruthTable(p.TTName)
	if err != nil {
		notify.Fatal(err)
	}

	// Begin a performance profile.
	if p.ProfName != "" {
		f, err := os.Create(p.ProfName)
		if err != nil {
			notify.Fatal(err)
		}
		defer f.Close()
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	// Solve for the QUBO coefficients that maximize the gap.
	tt, gap, vals, coeffs, ok := FindCoefficients(&p, tt)
	if !ok {
		notify.Fatal("Failed to solve for the QUBO coefficients")
	}
	fmt.Printf("QUBO = %s\n", coeffsToOctave(tt.NCols, coeffs))
	fmt.Printf("Gap = %v\n", gap)
	outputEvaluation(&p, tt, vals)
}
