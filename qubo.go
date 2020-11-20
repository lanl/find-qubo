// This file tries to find the QUBO parameters that implement a given
// truth table.

package main

import (
	"fmt"
	"math"
	"math/rand"
	"strings"
	"time"

	"github.com/MaxHalford/eaopt"
	"gonum.org/v1/gonum/mat"
)

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
	Gap    float64     // Difference in value between the maximum valid row and the minimum invalid row
}

// NewSPSOQUBO runs a quick particle-swarm optimization to choose initial
// QUBO coefficients for further genetic-algorithm optimization.
func NewSPSOQUBO(p *Parameters, rng *rand.Rand) *QUBO {
	spso, err := eaopt.NewDefaultSPSO()
	if err != nil {
		notify.Fatal(err)
	}
	spso.Min = math.Max(p.MinL, p.MinQ)
	spso.Max = math.Min(p.MaxL, p.MaxQ)
	cfs, _, err := spso.Minimize(func(cfs []float64) float64 {
		qubo := &QUBO{
			Params: p,
			Coeffs: cfs,
			Gap:    -math.MaxFloat64,
		}
		bad, _ := qubo.Evaluate()
		return bad
	}, uint((p.NCols*(p.NCols+1))/2))
	if err != nil {
		notify.Fatal(err)
	}
	return &QUBO{
		Params: p,
		Coeffs: cfs,
		Gap:    -math.MaxFloat64,
	}
}

// EvaluateAllInputs multiplies the QUBO by each input column in turn (i.e.,
// x'*Q*x for all x).
func (q *QUBO) EvaluateAllInputs() []float64 {
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

// SelectValidRows is a helper function for Evaluate that selects at most one
// row in each batch of 2^a (with a ancillae) to treat as valid.
func (q *QUBO) SelectValidRows(vals []float64) []bool {
	// Iterate over each base row (i.e., rows if there were no ancillae).
	p := q.Params
	nRows := len(p.TT)
	valids := make([]bool, nRows)
	bSize := 1 << p.NAnc // Number of rows to consider as a single batch
	for base := 0; base < nRows; base += bSize {
		// Select the smallest-valued row as the valid row.
		if !p.TT[base] {
			continue // No rows in the batch should be considered valid.
		}
		var vRow int // Single valid row
		minVal := math.MaxFloat64
		for ofs := 0; ofs < bSize; ofs++ {
			r := base + ofs
			if vals[r] < minVal {
				minVal = vals[r]
				vRow = r
			}
		}
		valids[vRow] = true
	}
	return valids
}

// Evaluate computes the badness of a set of coefficients.
func (q *QUBO) Evaluate() (float64, error) {
	// Find the minimum output across all inputs.
	p := q.Params
	vals := q.EvaluateAllInputs()
	minVal := math.MaxFloat64
	for _, v := range vals {
		minVal = math.Min(minVal, v)
	}

	// Find the maximum valid output and the minimum invalid output.
	maxValid := -math.MaxFloat64
	minInvalid := math.MaxFloat64
	isValid := q.SelectValidRows(vals)
	for r, v := range vals {
		if isValid[r] {
			maxValid = math.Max(maxValid, v)
		} else {
			minInvalid = math.Min(minInvalid, v)
		}
	}
	q.Gap = minInvalid - maxValid

	// Penalize valid rows in the truth table that produced a non-minimal
	// value and invalid rows that produced a value better than any valid
	// row.
	bad := 0.0
	wt := float64(len(p.TT)) * math.Max(p.MaxL, p.MaxQ)
	for r, v := range vals {
		switch {
		case isValid[r]:
			// Valid row: Penalize according to the value's amount
			// above the global minimal value.
			bad += math.Pow(v-minVal, 2.0)
		case v <= maxValid:
			// Invalid row with a value less than or equal to the
			// maximum valid value: severely penalize according to
			// the value's amount below the maximum valid value.
			bad += (math.Pow(v-maxValid, 2.0) + 1) * wt
		default:
			// Invalid row with a value greater than the maximum
			// valid value: No penalty.
		}
	}

	// Add another penalty for small gaps, but only after we've
	// successfully separated valid from invalid rows.
	bad = (bad + 1.0) * p.MaxGap
	if p.RewardGap && minInvalid > maxValid {
		bad -= minInvalid - maxValid // Reward large gaps.
	}
	return bad, nil
}

// mutateRandomize mutates a single coefficient at random.
func (q *QUBO) mutateRandomize(rng *rand.Rand) {
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

// mutateRandomizeAll mutates all coefficients at random.
func (q *QUBO) mutateRandomizeAll(rng *rand.Rand) {
	p := q.Params
	for c := range q.Coeffs {
		if c < p.NCols {
			// Linear coefficient
			q.Coeffs[c] = rng.Float64()*(p.MaxL-p.MinL) + p.MinL
		} else {
			// Quadratic coefficient
			q.Coeffs[c] = rng.Float64()*(p.MaxQ-p.MinQ) + p.MinQ
		}
	}
}

// mutateRound rounds all coefficients to the nearest N.
func (q *QUBO) mutateRound(rng *rand.Rand) {
	ns := [...]float64{0x1p-8, 0x1p-16, 0x1p-24}
	n := ns[rng.Intn(len(ns))]
	for i, c := range q.Coeffs {
		q.Coeffs[i] = math.Round(c/n) * n
	}
}

// mutateFlipSign negates a single coefficient at random.
func (q *QUBO) mutateFlipSign(rng *rand.Rand) {
	c := rng.Intn(len(q.Coeffs))
	q.Coeffs[c] = -q.Coeffs[c]
}

// mutateCopy copies one coefficient to another.
func (q *QUBO) mutateCopy(rng *rand.Rand) {
	// Select two coefficients at random.
	nc := len(q.Coeffs)
	c1 := rng.Intn(nc)
	c2 := (rng.Intn(nc-1) + c1 + 1) % nc // Different from c1

	// Sometimes replace these with the two nearest coefficients in
	// magnitude that are not already equal.
	if rng.Intn(5) == 0 {
		closest := math.MaxFloat64
		for i := 0; i < nc-1; i++ {
			v1 := math.Abs(q.Coeffs[i])
			for j := i + 1; j < nc; j++ {
				v2 := math.Abs(q.Coeffs[j])
				d := math.Abs(v1 - v2)
				if d != 0.0 && d < closest {
					closest = d
					c1 = i
					c2 = j
				}
			}
		}
	}
	q.Coeffs[c1] = math.Copysign(q.Coeffs[c2], q.Coeffs[c1]) // Copy only the magnitude.
}

// mutateNudge slightly modifies a single coefficient.
func (q *QUBO) mutateNudge(rng *rand.Rand) {
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
func (q *QUBO) Mutate(rng *rand.Rand) {
	r := rng.Intn(100)
	switch {
	case r == 0:
		// Randomize all coefficients (rare).
		q.mutateRandomizeAll(rng)
	case r == 1:
		// Round all coefficients (rare).
		q.mutateRound(rng)
	case r < 5:
		// Negate a coefficient (fairly rare).
		q.mutateFlipSign(rng)
	case r < 5+30:
		// Copy one coefficient to another.
		q.mutateCopy(rng)
	case r < 5+30+20:
		// Randomize a single coefficient.
		q.mutateRandomize(rng)
	default:
		// Slightly modify a single coefficient.
		q.mutateNudge(rng)
	}
}

// Crossover randomly blends the coefficients of two QUBOs.
func (q *QUBO) Crossover(g eaopt.Genome, rng *rand.Rand) {
	q2 := g.(*QUBO)
	eaopt.CrossUniformFloat64(q.Coeffs, q2.Coeffs, rng)
}

// Clone returns a copy of a QUBO.
func (q *QUBO) Clone() eaopt.Genome {
	cfs := make([]float64, len(q.Coeffs))
	copy(cfs, q.Coeffs)
	return &QUBO{
		Params: q.Params,
		Coeffs: cfs,
		Gap:    q.Gap,
	}
}

// Rescale scales all coefficients so the maximum is as large as possible.
func (q *QUBO) Rescale() {
	// Find the maximal linear term.
	p := q.Params
	nc := p.NCols
	maxLin := -math.MaxFloat64
	for i := 0; i < nc; i++ {
		maxLin = math.Max(maxLin, math.Abs(q.Coeffs[i]))
	}

	// Find the maximal quadratic term.
	maxQuad := -math.MaxFloat64
	for i := nc; i < len(q.Coeffs); i++ {
		maxQuad = math.Max(maxQuad, math.Abs(q.Coeffs[i]))
	}

	// Scale all coefficients equally.
	maxL := math.Min(p.MaxL, -p.MinL)
	maxQ := math.Min(p.MaxQ, -p.MinQ)
	scale := math.Min(maxL/maxLin, maxQ/maxQuad)
	for i := range q.Coeffs {
		q.Coeffs[i] *= scale
	}
}

// AsOctaveMatrix returns the coefficients as a string that can be pasted into
// GNU Octave or MATLAB.
func (q *QUBO) AsOctaveMatrix() string {
	p := q.Params
	i := p.NCols
	oct := make([]string, p.NCols)
	for r := 0; r < p.NCols; r++ {
		row := make([]string, p.NCols)
		for c := 0; c < p.NCols; c++ {
			switch {
			case c < r:
				// Not in upper triangle
				row[c] = "0"
			case c == r:
				// Linear term
				row[c] = fmt.Sprint(q.Coeffs[c])
			default:
				// Quadratic term
				row[c] = fmt.Sprint(q.Coeffs[i])
				i++
			}
		}
		oct[r] = strings.Join(row, " ")
	}
	return "[" + strings.Join(oct, " ; ") + "]"
}

// OptimizeCoeffs tries to find the coefficients that best represent the given
// truth table, the corresponding badness, and the number of generations
// evolved.  It aborts on error.
func OptimizeCoeffs(p *Parameters) (*QUBO, float64, uint) {
	// Create a genetic-algorithm object.
	cfg := eaopt.NewDefaultGAConfig()
	cfg.NGenerations = 1000000
	cfg.Model = eaopt.ModGenerational{
		Selector:  eaopt.SelElitism{},
		MutRate:   0.85,
		CrossRate: 0.50,
	}
	prevBest := math.MaxFloat64 // Least badness seen so far
	startTime := time.Now()     // Current time
	prevReport := startTime     // Last time we reported our status
	separatedGen := -1          // First generation at which an HOF QUBO separated valid from invalid rows
	cfg.Callback = func(ga *eaopt.GA) {
		hof := ga.HallOfFame[0]
		bad := hof.Fitness
		if bad < prevBest && time.Since(prevReport) > 3*time.Second {
			// Report when we have a new least badness but not more
			// than once per second.
			status.Printf("Least badness = %.10g after %d generations and %.1fs", bad, ga.Generations, ga.Age.Seconds())
			qubo := hof.Genome.(*QUBO)
			status.Printf("    Best coefficients = %v", qubo.Coeffs)
			status.Printf("    Matrix form = %s", qubo.AsOctaveMatrix())
			if qubo.Gap >= 0.0 {
				status.Printf("    Valid/invalid gap = %v", qubo.Gap)
			}
			if qubo.Gap >= 0.0 && separatedGen == -1 {
				status.Print("All valid rows finally have lower values than all invalid rows!")
				status.Printf("Running for %d more generations in attempt to increase the valid/invalid gap", p.GapIters)
				separatedGen = int(ga.Generations)
				p.RewardGap = true
			}
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
		return separatedGen >= 0 && int(ga.Generations)-separatedGen > p.GapIters
	}
	ga, err := cfg.NewGA()
	if err != nil {
		notify.Fatal(err)
	}

	// Run the genetic algorithm.
	err = ga.Minimize(func(rng *rand.Rand) eaopt.Genome {
		return NewSPSOQUBO(p, rng)
	})
	if err != nil {
		notify.Fatal(err)
	}

	// Return the best coefficients we found.
	hof := ga.HallOfFame[0]
	return hof.Genome.(*QUBO), hof.Fitness, ga.Generations
}
