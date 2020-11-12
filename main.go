/*
Find parameters for a QUBO given a truth table.
*/

package main

import (
	"log"
	"os"
)

// notify is used to output error messages.
var notify *log.Logger

func main() {
	// Initialize program parameters.
	notify = log.New(os.Stderr, os.Args[0]+": ", 0)
	var p Parameters
	ParseCommandLine(&p)

	// Read the input file.
	tt, nc := ReadTruthTable(&p)
	p.TT = tt
	p.NCols = nc

	// Precompute a matrix with all possible 0/1 columns.
	p.AllCols = AllPossibleColumns(p.NCols)

	// Temporary
	notify.Printf("TT = %v", tt)
	qubo := QUBO{
		Params: &p,
		Coeffs: []float64{0, 0, 3.0, 1.0, -2.0, -2.0},
	}
	bad, _ := qubo.Evaluate()
	notify.Printf("Bad = %v", bad)
}
