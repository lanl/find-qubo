// This file defines program parameters and routines for initializing them.

package main

import (
	"flag"
	"fmt"
	"math"
	"os"

	"gonum.org/v1/gonum/mat"
)

// Parameters is a collection of all program parameters.
type Parameters struct {
	TTName       string     // Name of the input truth-table file
	MinQ         float64    // Minimum quadratic coefficient
	MaxQ         float64    // Maximum quadratic coefficient
	MinL         float64    // Minimum linear coefficient
	MaxL         float64    // Maximum linear coefficient
	TT           TruthTable // The truth-table proper
	NCols        int        // Number of columns in the truth table, including ancillae
	NAnc         int        // Number of ancilla columns
	AllCols      *mat.Dense // Matrix with all 2^n columns for n rows
	Balance      bool       // true=also consider valid-row balance; false=consider only gap
	ZeroGen      int        // Iteration at which the gapped dropped to near-zero without going below zero
	ZeroGenLen   int        // Number of generations to allow the GA to have a near-zero gap
	SeparatedGen int        // Iteration at which invalid and valid rows finally separated
	MaxGap       float64    // Loose upper bound on the maximum gap
	GapIters     int        // Number of iterations to perform to increase the valid/invalid gap
	RewardGap    bool       // true=reward large valid/invalid gaps (later iterations); false=ignore them (early iterations)
}

// PrepareGAParameters initializes parameters needed for the next run of the
// genetic algorithm.
func PrepareGAParameters(p *Parameters) {
	// We haven't yet separated valid from invalid rows or encountered
	// all-zero coefficients.
	p.SeparatedGen = -1
	p.ZeroGen = -1

	// Read the input file.
	tt, nc := ReadTruthTable(p)
	p.TT = tt
	p.NCols = nc

	// Precompute a matrix with all possible 0/1 columns.
	p.AllCols = AllPossibleColumns(p.NCols)

	// Precompute the maximum gap, rounded up to a power of 10.
	p.MaxGap = findMaxGap(p)
}

// findMaxGap computes a loose, worst-case estimate of the valid/invalid gap.
func findMaxGap(p *Parameters) float64 {
	nc := float64((p.NCols * (p.NCols + 1)) / 2) // Number of coefficients
	maxGap := math.Max(float64(p.MaxL-p.MinL), float64(p.MaxQ-p.MinQ)) * nc
	maxGap = math.Pow(10.0, math.Ceil(math.Log10(maxGap))) // Round up to a power of 10
	return maxGap
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
	flag.BoolVar(&p.Balance, "balance", false, "Try harder to balance the values of valid rows")
	flag.IntVar(&p.NAnc, "ancillae", 0, "Number of ancilla columns to add")
	flag.IntVar(&p.ZeroGenLen, "zero-gen-len", 5000, "Number of generations to allow the GA to remain stuck with a small but negative invalid-valid gap")
	flag.IntVar(&p.GapIters, "gap-iters", 100000, "Number of extra iterations to perform to increase the invalid-valid gap")
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
	case p.GapIters < 0:
		notify.Fatal("--gap-iters must be non-negative")
	case p.ZeroGenLen < 0:
		notify.Fatal("--zero-gen-len must be non-negative")
	}
}
