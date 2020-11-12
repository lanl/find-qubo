// This file defines program parameters and routines for initializing
// them from the command line.

package main

import (
	"flag"
	"fmt"
	"os"

	"gonum.org/v1/gonum/mat"
)

// Parameters is a collection of all program parameters.
type Parameters struct {
	TTName  string     // Name of the input truth-table file
	MinJ    float64    // Minimum quadratic coefficient
	MaxJ    float64    // Maximum quadratic coefficient
	MinH    float64    // Minimum linear coefficient
	MaxH    float64    // Maximum linear coefficient
	TT      TruthTable // The truth-table proper
	NCols   int        // Number of columns in the truth table
	AllCols *mat.Dense // Matrix with all 2^n columns for n rows
}

// ParseCommandLine parses parameters from the command line.
func ParseCommandLine(p *Parameters) {
	// Parse the command line.
	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "Usage: %s [<options>] [<input.tt>]\n", os.Args[0])
		flag.PrintDefaults()
	}
	flag.Float64Var(&p.MinJ, "jmin", -1.0, "Minimum quadratic coefficient")
	flag.Float64Var(&p.MaxJ, "jmax", 1.0, "Maximum quadratic coefficient")
	flag.Float64Var(&p.MinH, "hmin", -2.0, "Minimum linear coefficient")
	flag.Float64Var(&p.MaxH, "hmax", 2.0, "Maximum linear coefficient")
	flag.Parse()
	if flag.NArg() >= 1 {
		p.TTName = flag.Arg(0)
	}

	// Validate the arguments.
	if p.MinJ >= p.MaxJ {
		notify.Fatal("--jmin must specify a value that is less than -jmax")
	}
	if p.MinH >= p.MaxH {
		notify.Fatal("--hmin must specify a value that is less than -hmax")
	}
}
