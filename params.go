// This file defines program parameters and routines for initializing them.

package main

import (
	"flag"
	"fmt"
	"os"
)

// Parameters is a collection of all program parameters.
type Parameters struct {
	TTName      string  // Name of the input truth-table file
	MinQ        float64 // Minimum quadratic coefficient
	MaxQ        float64 // Maximum quadratic coefficient
	MinL        float64 // Minimum linear coefficient
	MaxL        float64 // Maximum linear coefficient
	RoundTo     float64 // Value to which to round all coefficients
	MaxAncillae uint    // Maximum number of additional variables we're allowed to add
	ProfName    string  // Name of a pprof performance-profile file
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
	flag.Float64Var(&p.RoundTo, "round", 0, "Value to which to round coefficients or 0 for no rounding")
	flag.UintVar(&p.MaxAncillae, "max-ancillae", 10, "Maximum number of ancilllary variables the program is allowed to add")
	flag.StringVar(&p.ProfName, "profile", "", "Name of a pprof performance file to write")
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
	case p.RoundTo < 0.0:
		notify.Fatal("--round must be non-negative")
	}
}
