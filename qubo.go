// This file tries to find the QUBO parameters that implement a given
// truth table.

package main

import (
	"math"
	"math/rand"
	"time"

	"github.com/MaxHalford/eaopt"
	"gonum.org/v1/gonum/mat"
)

// GoodEnoughBadness is a badness level we consider good enough to terminate
// optimization.
const GoodEnoughBadness = 1e-10

// AllPossibleColumns returns a matrix containing all possible columns
// for a given number of rows that contain only 0s and 1s.
func AllPossibleColumns(nr int) *mat.Dense {
	nc := 1 << nr
	all := mat.NewDense(nr, nc, make([]float64, nr*nc))
	for r := 0; r < nr; r++ {
		for c := 0; c < nc; c++ {
			if c&(1<<r) != 0 {
				all.Set(nr-r-1, c, 1.0)
			}
		}
	}
	return all
}

// A QUBO represents a QUBO matrix whose values we're solving for.
type QUBO struct {
	Params *Parameters // Pointer to global program parameters
	Coeffs []float64   // List of linear followed by quadratic coefficients
}

// NewRandomQUBO returns a QUBO with random coefficients.
func NewRandomQUBO(p *Parameters, rng *rand.Rand) QUBO {
	nCfs := p.NCols * (p.NCols + 1) / 2
	cfs := make([]float64, nCfs)
	for c := range cfs {
		if c < p.NCols {
			// Linear coefficient
			cfs[c] = rng.Float64()*(p.MaxL-p.MinL) + p.MinL
		} else {
			// Quadratic coefficient
			cfs[c] = rng.Float64()*(p.MaxQ-p.MinQ) + p.MinQ
		}
	}
	return QUBO{
		Params: p,
		Coeffs: cfs,
	}
}

// EvaluateAllInputs multiplies the QUBO by each input column in turn (i.e.,
// x'*Q*x for all x).
func (q QUBO) EvaluateAllInputs() []float64 {
	// Convert the coefficients to an upper-triangular matrix.
	p := q.Params // Global parameters
	n := p.NCols  // Number of columns and rows
	Q := mat.NewSymDense(n, make([]float64, n*n))
	for i := 0; i < n; i++ {
		Q.SetSym(i, i, q.Coeffs[i]) // Put linear coefficients on the diagonal.
	}
	k := n
	for i := 0; i < n-1; i++ {
		for j := i + 1; j < n; j++ {
			Q.SetSym(i, j, q.Coeffs[k]/2.0) // Put quadratic coefficients off the diagonal.
			k++
		}
	}

	// Evaluate the matrix on all 2^n possible inputs.
	vals := make([]float64, len(p.TT))
	for r := range p.TT {
		col := p.AllCols.ColView(r)
		m := mat.NewDense(1, 1, []float64{0.0})
		m.Product(col.T(), Q, col)
		vals[r] = m.At(0, 0)
	}
	return vals
}

// Evaluate computes the badness of a set of coefficients.
func (q QUBO) Evaluate() (float64, error) {
	// Find the minimum output across all inputs.
	p := q.Params
	vals := q.EvaluateAllInputs()
	minVal := math.MaxFloat64
	for _, v := range vals {
		minVal = math.Min(minVal, v)
	}

	// Penalize valid rows in the truth table that produced a
	// non-minimal value.
	bad := 0.0
	for r, v := range vals {
		if !p.TT[r] {
			continue
		}
		bad += math.Pow(v-minVal, 2.0)
	}
	return bad, nil
}

// mutateRandomize mutates a single coefficient at random.
func (q QUBO) mutateRandomize(rng *rand.Rand) {
	p := q.Params
	c := rng.Intn(len(q.Coeffs))
	if c < p.NCols {
		// Linear coefficient
		q.Coeffs[c] = rng.Float64()*(p.MaxL-p.MinL) + p.MinL
	} else {
		// Quadratic coefficient
		q.Coeffs[c] = rng.Float64()*(p.MaxQ-p.MinQ) + p.MinQ
	}
}

// mutateRound rounds all coefficients to the nearest N.
func (q QUBO) mutateRound(n float64) {
	for i, c := range q.Coeffs {
		q.Coeffs[i] = math.Round(c/n) * n
	}
}

// mutateCopy copies one coefficient to another.
func (q QUBO) mutateCopy(rng *rand.Rand) {
	nc := len(q.Coeffs)
	c1 := rng.Intn(nc)
	c2 := (rng.Intn(nc-1) + c1 + 1) % nc                     // Different from c1
	q.Coeffs[c1] = math.Copysign(q.Coeffs[c2], q.Coeffs[c1]) // Copy only the magnitude.
}

// mutateNudge slightly modifies a single coefficient.
func (q QUBO) mutateNudge(rng *rand.Rand) {
	const nudgeAmt = 1e-4
	p := q.Params
	c := rng.Intn(len(q.Coeffs))
	cf := q.Coeffs[c] + rng.Float64()*nudgeAmt*2 - nudgeAmt
	if c < p.NCols {
		// Linear coefficient
		cf = math.Min(math.Max(cf, p.MinL), p.MaxL)
	} else {
		// Quadratic coefficient
		cf = math.Min(math.Max(cf, p.MinQ), p.MaxQ)
	}
	q.Coeffs[c] = cf
}

// Mutate mutates the QUBO's coefficients.
func (q QUBO) Mutate(rng *rand.Rand) {
	r := rng.Intn(100)
	switch {
	case r == 0:
		// Round all coefficients (rare).
		q.mutateRound(1e-10)
	case r < 20:
		// Copy one coefficient to another.
		q.mutateCopy(rng)
	case r < 20+30:
		// Randomize a single coefficient.
		q.mutateRandomize(rng)
	default:
		// Slightly modify a single coefficient.
		q.mutateNudge(rng)
	}
}

// Crossover randomly blends the coefficients of two QUBOs.
func (q QUBO) Crossover(g eaopt.Genome, rng *rand.Rand) {
	q2 := g.(QUBO)
	eaopt.CrossGNXFloat64(q.Coeffs, q2.Coeffs, 2, rng)
}

// Clone returns a copy of a QUBO.
func (q QUBO) Clone() eaopt.Genome {
	cfs := make([]float64, len(q.Coeffs))
	copy(cfs, q.Coeffs)
	return QUBO{
		Params: q.Params,
		Coeffs: cfs,
	}
}

// OptimizeCoeffs tries to find the coefficients that best represent the given
// truth table and the corresponding badness.  It aborts on error.
func OptimizeCoeffs(p *Parameters) (QUBO, float64) {
	// Create a genetic-algorithm object.
	cfg := eaopt.NewDefaultGAConfig()
	cfg.NGenerations = 1000000
	cfg.Model = eaopt.ModGenerational{
		Selector:  eaopt.SelElitism{},
		MutRate:   0.75,
		CrossRate: 0.75,
	}
	prevBest := math.MaxFloat64 // Least badness seen so far
	startTime := time.Now()     // Current time
	prevReport := startTime     // Last time we reported our status
	cfg.Callback = func(ga *eaopt.GA) {
		hof := ga.HallOfFame[0]
		bad := hof.Fitness
		if bad < prevBest {
			// Report when we have a new least badness.
			status.Printf("Least badness = %.10g after %d generations and %.1fs", bad, ga.Generations, ga.Age.Seconds())
			status.Printf("Best coefficients = %v", hof.Genome.(QUBO).Coeffs)
			prevBest = bad
			prevReport = time.Now()
			return
		}
		if time.Since(prevReport) > 5*time.Second {
			status.Printf("Working on generation %d at time %.1fs", ga.Generations, time.Since(startTime).Seconds())
			prevReport = time.Now()
		}
	}
	cfg.EarlyStop = func(ga *eaopt.GA) bool {
		return ga.HallOfFame[0].Fitness <= GoodEnoughBadness
	}
	ga, err := cfg.NewGA()
	if err != nil {
		notify.Fatal(err)
	}

	// Run the genetic algorithm.
	err = ga.Minimize(func(rng *rand.Rand) eaopt.Genome {
		return NewRandomQUBO(p, rng)
	})
	if err != nil {
		notify.Fatal(err)
	}

	// Return the best coefficients we found.
	hof := ga.HallOfFame[0]
	return hof.Genome.(QUBO), hof.Fitness
}
