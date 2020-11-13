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

// status is used to output status messages.
var status *log.Logger

func main() {
	// Initialize program parameters.
	notify = log.New(os.Stderr, os.Args[0]+": ", 0)
	status = log.New(os.Stderr, "INFO: ", 0)
	var p Parameters
	ParseCommandLine(&p)

	// Read the input file.
	tt, nc := ReadTruthTable(&p)
	p.TT = tt
	p.NCols = nc

	// Precompute a matrix with all possible 0/1 columns.
	p.AllCols = AllPossibleColumns(p.NCols)

	// Try to find coefficients that represent the truth table.
	qubo, bad := OptimizeCoeffs(&p)
	status.Printf("Best coefficients = %v", qubo.Coeffs)
	status.Printf("Badness = %v", bad)
}
