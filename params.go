// This file defines program parameters and routines for initializing them.

package main

import (
	"flag"
	"fmt"
	"os"
)

// Parameters is a collection of all program parameters.
type Parameters struct {
	TTName      string            // Name of the input truth-table file
	MinQ        float64           // Minimum quadratic coefficient
	MaxQ        float64           // Maximum quadratic coefficient
	MinL        float64           // Minimum linear coefficient
	MaxL        float64           // Maximum linear coefficient
	RoundTo     float64           // Value to which to round all coefficients
	MaxAncillae uint              // Maximum number of additional variables we're allowed to add
	ProfName    string            // Name of a pprof performance-profile file
	Tolerance   float64           // Smallest-in-magnitude values for the LP solver to consider nonzero
	NumLPSolves uint64            // Tally of the number of LP solver invocations
	Approach    ReductionApproach // How to reduce the exponential search space
	Rank        int               // The current process's rank in the parallel computation
	NumRanks    int               // The total number of ranks in the parallel computation
}

// A ReductionApproach defines an approach to reduce the search space.
type ReductionApproach int

// These are the acceptable values for a ReductionApproach.
const (
	ReduceHeuristic     ReductionApproach = iota // Use a heuristic approach.
	ReduceBruteForce                             // Try all possibilities until one succeeds.
	ReduceBruteForceAll                          // Try all possibilities and tally successes/failures.
)

// String returns a ReductionApproach as a string.
func (ra *ReductionApproach) String() string {
	switch *ra {
	case ReduceHeuristic:
		return "heuristic"
	case ReduceBruteForce:
		return "brute-force"
	case ReduceBruteForceAll:
		return "full-brute-force"
	default:
		panic(fmt.Sprintf("unexpected ReductionApproach %d", *ra))
	}
}

// Set assigns a ReductionApproach from a string.
func (ra *ReductionApproach) Set(s string) error {
	switch s {
	case "heuristic":
		*ra = ReduceHeuristic
	case "brute-force":
		*ra = ReduceBruteForce
	case "full-brute-force":
		*ra = ReduceBruteForceAll
	default:
		return fmt.Errorf("unexpected reduction approach %q", s)
	}
	return nil
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
	flag.Float64Var(&p.Tolerance, "tolerance", 1e-10, "Smallest-in-magnitude values for the LP solver to consider nonzero")
	flag.Var(&p.Approach, "approach", `Approach to reducing the search space, one of "heuristic", "brute-force", or "full-brute-force" (default: heuristic)`)
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
	case p.Tolerance < 0.0:
		notify.Fatal("--tolerance must be non-negative")
	}
}
