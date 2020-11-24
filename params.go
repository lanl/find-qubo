// This file defines program parameters and routines for initializing them.

package main

import (
	"flag"
	"fmt"
	"os"

	"gonum.org/v1/gonum/mat"
)

// Parameters is a collection of all program parameters.
type Parameters struct {
	TTName     string     // Name of the input truth-table file
	MinQ       float64    // Minimum quadratic coefficient
	MaxQ       float64    // Maximum quadratic coefficient
	MinL       float64    // Minimum linear coefficient
	MaxL       float64    // Maximum linear coefficient
	TT         TruthTable // The truth-table proper
	NCols      int        // Number of columns in the truth table, including ancillae
	NAnc       int        // Number of ancilla columns
	AllCols    *mat.Dense // Matrix with all 2^n columns for n rows
	ZeroGen    int        // Iteration at which the gap appeared to reach an asymptote just below zero
	ZeroGenLen int        // Number of generations to allow the GA to have a near-zero negative gap
}

// PrepareGAParameters initializes parameters needed for the next run of the
// genetic algorithm.
func PrepareGAParameters(p *Parameters) {
	// We haven't yet encountered a near-zero negative gap.
	p.ZeroGen = -1

	// Precompute a matrix with all possible 0/1 columns.
	p.AllCols = AllPossibleColumns(p.NCols)
}

// ParseCommandLine parses parameters from the command line.
func ParseCommandLine(p *Parameters) {
	// Parse the command line.
	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "Usage: %s [<options>] [<input.tt>]\n", os.Args[0])
		flag.PrintDefaults()
	}
	flag.Float64Var(&p.MinQ, "qmin", -1.0, "Minimum quadratic coefficient")
	flag.Float64Var(&p.MaxQ, "qmax", 1.0, "Maximum quadratic coefficient")
	flag.Float64Var(&p.MinL, "lmin", -1.0, "Minimum linear coefficient")
	flag.Float64Var(&p.MaxL, "lmax", 1.0, "Maximum linear coefficient")
	flag.IntVar(&p.NAnc, "ancillae", 0, "Initial number of ancilla columns to add")
	flag.IntVar(&p.ZeroGenLen, "zero-gen-len", 5000, "Number of generations to allow the GA to remain stuck with a small but negative invalid-valid gap")
	flag.Parse()
	if flag.NArg() >= 1 {
		p.TTName = flag.Arg(0)
	}

	// Validate the arguments.
	switch {
	case p.MinQ >= p.MaxQ:
		notify.Fatal("--qmin must specify a value that is less than --qmax")
	case p.MinL >= p.MaxL:
		notify.Fatal("--lmin must specify a value that is less than --lmax")
	case p.MinL >= 0.0:
		notify.Fatal("--lmin must be negative")
	case p.MinQ >= 0.0:
		notify.Fatal("--qmin must be negative")
	case p.MaxL <= 0.0:
		notify.Fatal("--lmax must be positive")
	case p.MaxQ <= 0.0:
		notify.Fatal("--qmax must be positive")
	case p.NAnc < 0:
		notify.Fatal("--ancillae must be non-negative")
	case p.ZeroGenLen < 0:
		notify.Fatal("--zero-gen-len must be non-negative")
	}
}
