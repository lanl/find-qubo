// This file tries to find the QUBO parameters that implement a given
// truth table.

package main

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"strings"
	"time"

	"github.com/MaxHalford/eaopt"
	"github.com/lanl/clp"
	"gonum.org/v1/gonum/mat"
)

// NearZero defines a value close enough to zero to be considered a zero gap.
const NearZero = 1e-5

// A QUBO represents a QUBO matrix whose values we're solving for.
type QUBO struct {
	Params  *Parameters    // Pointer to global program parameters
	Coeffs  []float64      // List of linear followed by quadratic coefficients
	Gap     float64        // Difference in value between the maximum valid row and the minimum invalid row
	History map[string]int // Tally of genome modifications made
}

// NewSPSOQUBO runs a quick particle-swarm optimization to choose initial
// QUBO coefficients for further genetic-algorithm optimization.
func NewSPSOQUBO(p *Parameters, rng *rand.Rand) *QUBO {
	qubo := &QUBO{
		Params:  p,
		Coeffs:  coeffsSPSO(p, rng),
		Gap:     -math.MaxFloat64,
		History: make(map[string]int, 100),
	}
	return qubo
}

// NewRandomQUBO chooses initial QUBO coefficients at random.
func NewRandomQUBO(p *Parameters, rng *rand.Rand) *QUBO {
	qubo := &QUBO{
		Params:  p,
		Coeffs:  coeffsRandom(p, rng, 0.0),
		Gap:     -math.MaxFloat64,
		History: make(map[string]int, 100),
	}
	return qubo
}

// coeffsSPSO returns coefficients found from particle-swarm optimization.
// Empirical results indicate that these are likely to be a decent local
// minimum.  The function aborts on error.
func coeffsSPSO(p *Parameters, rng *rand.Rand) []float64 {
	spso, err := eaopt.NewDefaultSPSO()
	if err != nil {
		notify.Fatal(err)
	}
	spso.Min = math.Max(p.MinL, p.MinQ)
	spso.Max = math.Min(p.MaxL, p.MaxQ)
	bestCfs, _, err := spso.Minimize(func(cfs []float64) float64 {
		if math.IsNaN(cfs[0]) {
			// I don't know how we could have gotten here, but I
			// have seen an array of NaNs passed in.
			return math.MaxFloat64
		}
		qubo := &QUBO{
			Params:  p,
			Coeffs:  cfs,
			Gap:     -math.MaxFloat64,
			History: make(map[string]int, 100),
		}
		bad, _ := qubo.Evaluate()
		return bad
	}, uint((p.NCols*(p.NCols+1))/2))
	if err != nil {
		notify.Fatal(err)
	}
	return bestCfs
}

// coeffsRandom returns a completely random set of coefficients, optionally
// rounding them to a given value.
func coeffsRandom(p *Parameters, rng *rand.Rand, rt float64) []float64 {
	// Generate random coefficients.
	nc := p.NCols
	cfs := make([]float64, (nc*(nc+1))/2)
	for c := range cfs {
		if c < nc {
			// Linear coefficient
			cfs[c] = rng.Float64()*(p.MaxL-p.MinL) + p.MinL
		} else {
			// Quadratic coefficient
			cfs[c] = rng.Float64()*(p.MaxQ-p.MinQ) + p.MinQ
		}
	}

	// Optionally round each of them towards zero.
	if rt <= 0.0 {
		return cfs
	}
	for i, v := range cfs {
		if v >= 0.0 {
			cfs[i] = math.Floor(v/rt) * rt
		} else {
			cfs[i] = math.Ceil(v/rt) * rt
		}
	}
	return cfs
}

// coeffsBiased returns a set of coefficients that favors a single row of the
// truth table, ignoring all others.
func coeffsBiased(p *Parameters, rng *rand.Rand) []float64 {
	// Select a row at random.
	row := rng.Intn(len(p.TT))
	valid := p.TT[row]

	// Assign linear coefficients.
	nCfs := p.NCols * (p.NCols + 1) / 2
	cfs := make([]float64, nCfs)
	for b := 0; b < p.NCols; b++ {
		rb := p.NCols - b - 1 // Reverse bit order.
		v := (row >> rb) & 1
		switch {
		case valid && v == 0:
			cfs[b] = 1.0
		case valid && v == 1:
			cfs[b] = -1.0
		case !valid && v == 0:
			cfs[b] = -1.0
		case !valid && v == 1:
			cfs[b] = 1.0
		}
	}

	// Assign quadraric coefficients.
	i := p.NCols
	for b0 := 0; b0 < p.NCols-1; b0++ {
		rb0 := p.NCols - b0 - 1 // Reverse bit order.
		v0 := (row >> rb0) & 1
		for b1 := b0 + 1; b1 < p.NCols; b1++ {
			rb1 := p.NCols - b1 - 1 // Reverse bit order.
			v1 := (row >> rb1) & 1
			switch {
			case valid && v0 == v1:
				cfs[b0] += 1.0
				cfs[b1] += 1.0
				cfs[i] -= 2.0
			case valid && v0 != v1:
				cfs[b0] -= 1.0
				cfs[b1] -= 1.0
				cfs[i] += 2.0
			case !valid && v0 == v1:
				cfs[b0] -= 1.0
				cfs[b1] -= 1.0
				cfs[i] += 2.0
			case !valid && v0 != v1:
				cfs[b0] += 1.0
				cfs[b1] += 1.0
				cfs[i] -= 2.0
			}
			i++
		}
	}
	return cfs
}

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

	// Penalize misordered rows.
	bad := 0.0
	for r, v := range vals {
		const epsilon = 1e-5 // Small amount to add to discourage all-zero coefficients
		switch {
		case isValid[r] && v >= minInvalid:
			bad += math.Pow(epsilon+v-minInvalid, 2.0)
		case !isValid[r] && v <= maxValid:
			bad += math.Pow(epsilon+v-maxValid, 2.0)
		}
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

// mutateReplaceAll mutates all coefficients at random.
func (q *QUBO) mutateReplaceAll(rng *rand.Rand) {
	// Select a set of coefficients to consider.
	var cfs []float64
	switch rng.Intn(4) {
	case 0:
		// Local minimum found by particle-swarm optimization
		cfs = coeffsSPSO(q.Params, rng)
	case 1:
		// Coefficients biased to favor a single row
		cfs = coeffsBiased(q.Params, rng)
	case 2:
		// Coefficients biased to favor two rows
		cf1 := coeffsBiased(q.Params, rng)
		cf2 := coeffsBiased(q.Params, rng)
		for i, v := range cf2 {
			cf1[i] += v
		}
		cfs = cf1
	case 3:
		// Completely random coefficients, either rounded or not
		rtChoices := [...]float64{0.0, 0.03125, 0.0625, 0.125, 0.25, 0.5}
		rt := rtChoices[rng.Intn(len(rtChoices))]
		cfs = coeffsRandom(q.Params, rng, rt)
	default:
		panic("Unexpected option in mutateReplaceAll")
	}
	q.Coeffs = cfs
}

// mutateFlipSign negates a single coefficient at random.
func (q *QUBO) mutateFlipSign(rng *rand.Rand) {
	c := rng.Intn(len(q.Coeffs))
	q.Coeffs[c] = -q.Coeffs[c]
}

// mutateCopy copies one coefficient to another, optionally doubling one of the
// values.
func (q *QUBO) mutateCopy(rng *rand.Rand) {
	// Select two coefficients.
	nc := len(q.Coeffs)
	var c1, c2 int
	switch rng.Intn(5) {
	case 0:
		// Select the two nearest coefficients in magnitude that are
		// not already equal.
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
	default:
		// Select two coefficients at random.
		c1 = rng.Intn(nc)
		c2 = (rng.Intn(nc-1) + c1 + 1) % nc // Different from c1
	}

	// Copy one coefficient to the other either as is or halved in
	// magnitude.
	q.Coeffs[c1] = math.Copysign(q.Coeffs[c2], q.Coeffs[c1])
	if rng.Intn(3) == 0 {
		q.Coeffs[c1] /= 2.0
	}
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

// mutateReplaceLargest randomly replaces the largest-in-magnitude coefficient.
func (q *QUBO) mutateReplaceLargest(rng *rand.Rand) {
	// Find the largest coefficient.
	big := 0.0 // Largest magnitude coefficient (absolute value)
	idx := 0   // Index of the above
	for i, cf := range q.Coeffs {
		cf = math.Abs(cf)
		if cf > big {
			big = cf
			idx = i
		}
	}

	// Replace the coefficient with a random value, preserving its sign.
	p := q.Params
	var r float64
	if idx < p.NCols {
		// Linear coefficient
		r = rng.Float64()*(p.MaxL+p.MinL) - p.MinL
	} else {
		// Quadratic coefficient
		r = rng.Float64()*(p.MaxQ+p.MinQ) - p.MinQ
	}
	q.Coeffs[idx] = math.Copysign(r, q.Coeffs[idx])
}

// Mutate mutates the QUBO's coefficients.
func (q *QUBO) Mutate(rng *rand.Rand) {
	r := rng.Intn(100)
	switch {
	case r == 0:
		// Replace all coefficients (rare).
		q.mutateReplaceAll(rng)
		q.History["mutateReplaceAll"]++
	case r < 5:
		// Negate a coefficient (fairly rare).
		q.mutateFlipSign(rng)
		q.History["mutateFlipSign"]++
	case r < 5+20:
		// Copy one coefficient to another.
		q.mutateCopy(rng)
		q.History["mutateCopy"]++
	case r < 5+20+20:
		// Randomize a single coefficient.
		q.mutateRandomize(rng)
		q.History["mutateRandomize"]++
	case r < 5+20+20+20:
		// Replace the largest coefficient.
		q.mutateReplaceLargest(rng)
		q.History["mutateReplaceLargest"]++
	default:
		// Slightly modify a single coefficient.
		q.mutateNudge(rng)
		q.History["mutateNudge"]++
	}
}

// Crossover randomly blends the coefficients of two QUBOs.
func (q *QUBO) Crossover(g eaopt.Genome, rng *rand.Rand) {
	// Merge the two histories.  We use max instead of division to avoid
	// roughly doubling the values each iteration.
	q2 := g.(*QUBO)
	h1 := make(map[string]int, len(q.History)+len(q2.History))
	h2 := make(map[string]int, len(q.History)+len(q2.History))
	for k, v := range q.History {
		h1[k] = v
	}
	for k, v := range q2.History {
		if v > h1[k] {
			h1[k] = v
		}
	}
	for k, v := range h1 {
		h2[k] = v
	}
	q.History = h1
	q2.History = h2

	// Cross over the coefficients.
	eaopt.CrossUniformFloat64(q.Coeffs, q2.Coeffs, rng)
}

// Clone returns a copy of a QUBO.
func (q *QUBO) Clone() eaopt.Genome {
	// Copy the coefficients.
	cfs := make([]float64, len(q.Coeffs))
	copy(cfs, q.Coeffs)

	// Copy the history.
	hist := make(map[string]int, len(q.History))
	for k, v := range q.History {
		hist[k] = v
	}

	// Create and return a copy of the QUBO.
	return &QUBO{
		Params:  q.Params,
		Coeffs:  cfs,
		Gap:     q.Gap,
		History: hist,
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
	var scale float64
	switch {
	case maxLin == 0.0 && maxQuad == 0.0:
		// Both maxima are zero.
	case maxLin == 0.0:
		// All linear coefficients are zero.
		scale = maxQ / maxQuad
	case maxQuad == 0.0:
		// All quadratic coefficients are zero.
		scale = maxL / maxLin
	default:
		// Common case: Neither maximum is zero.
		scale = math.Min(maxL/maxLin, maxQ/maxQuad)
	}
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

// MakeGACallback returns a callback function to be executed after each
// generation of the genetic algorithm.
func MakeGACallback(p *Parameters) func(ga *eaopt.GA) {
	prevBest := math.MaxFloat64 // Least badness seen so far
	startTime := time.Now()     // Current time
	prevReport := startTime     // Last time we reported our status
	return func(ga *eaopt.GA) {
		hof := ga.HallOfFame[0]
		bad := hof.Fitness
		if bad < prevBest && time.Since(prevReport) > 3*time.Second {
			// Report when we have a new least badness but not more
			// than once every 3 seconds.
			qubo := hof.Genome.(*QUBO)
			status.Printf("Badness = %.10g (gap = %.1e) after %d generations and %.1fs", bad, qubo.Gap, ga.Generations, ga.Age.Seconds())

			// If the gap has been near zero for a long time
			// without crossing above zero, the problem is likely
			// unsolvable given the current number of ancillary
			// variables.
			switch {
			case p.ZeroGen == -1 && qubo.Gap <= 0.0 && -qubo.Gap < NearZero:
				// Small, negative gap: start the death clock
				// ticking.
				p.ZeroGen = int(ga.Generations)
			case p.ZeroGen != -1 && qubo.Gap > 0.0:
				// Any positive gap: stop the death clock.
				p.ZeroGen = -1
			}

			// Record when we last reported our status.
			prevBest = bad
			prevReport = time.Now()
			return
		}

		// In the case of no progress, provide a heartbeat message
		// every 5 seconds so the user knows we're still working.
		if time.Since(prevReport) > 5*time.Second {
			status.Printf("Working on generation %d at time %.1fs", ga.Generations, time.Since(startTime).Seconds())
			prevReport = time.Now()
		}
	}
}

// OptimizeCoeffs tries to find the coefficients that best represent the given
// truth table, the corresponding badness, and the number of generations
// evolved.  It aborts on error.
func OptimizeCoeffs(p *Parameters) (*QUBO, float64, uint) {
	// Create a genetic-algorithm object.
	cfg := eaopt.NewDefaultGAConfig()
	cfg.NGenerations = 1000000
	cfg.Model = eaopt.ModGenerational{
		Selector:  eaopt.SelTournament{NContestants: 3},
		MutRate:   0.85,
		CrossRate: 0.50,
	}
	cfg.PopSize = 30
	cfg.NPops = uint(runtime.NumCPU())
	cfg.Migrator = eaopt.MigRing{NMigrants: 5}
	cfg.MigFrequency = 10000
	cfg.Callback = MakeGACallback(p)
	cfg.EarlyStop = func(ga *eaopt.GA) bool {
		qubo := ga.HallOfFame[0].Genome.(*QUBO)
		if qubo.Gap > 0.0 {
			// We found a valid solution, meaning that all valid
			// rows evaluate to a smaller number than all invalid
			// rows.
			return true
		}
		if p.ZeroGen >= 0 && int(ga.Generations)-p.ZeroGen > p.ZeroGenLen {
			// A solution appears to be impossible given the current
			// number of ancillae.
			return true
		}
		return false // Keep running.
	}
	ga, err := cfg.NewGA()
	if err != nil {
		notify.Fatal(err)
	}

	// Run the genetic algorithm.
	err = ga.Minimize(func(rng *rand.Rand) eaopt.Genome {
		// Usually generate coefficients at random, but sometimes get
		// them from particle-swarm optimization.
		if rng.Intn(10) == 0 {
			return NewSPSOQUBO(p, rng)
		}
		return NewRandomQUBO(p, rng)
	})
	if err != nil {
		notify.Fatal(err)
	}

	// Return the best coefficients we found.
	hof := ga.HallOfFame[0]
	return hof.Genome.(*QUBO), hof.Fitness, ga.Generations
}

// LPReoptimize use a linear-programming algorithm to re-optimize a QUBO's
// coefficients in search of perfect balance of all valid rows.  This method
// returns a success code.
func (q *QUBO) LPReoptimize(isValid []bool) bool {
	// Define multipliers for each linear term's coefficients.  Each
	// multiplier will be either 0.0 (omitted) or 1.0.
	p := q.Params
	mat := clp.NewPackedMatrix()
	nc := len(q.Coeffs) // Number of coefficients
	nr := len(p.TT)     // Number of rows
	for c := 0; c < p.NCols; c++ {
		// Append one multiplier per truth-table row.
		cc := p.NCols - c - 1
		mults := make([]clp.Nonzero, 0, nr/2)
		for r := 0; r < nr; r++ {
			if (r>>cc)&1 == 1 {
				mults = append(mults, clp.Nonzero{
					Index: r,
					Value: 1.0,
				})
			}
		}

		// Append the complete column to the matrix.
		mat.AppendColumn(mults)
	}

	// Define multipliers for each quadratic term's coefficients.  Each
	// multiplier will be either 0.0 (omitted) or 1.0.
	for c1 := 0; c1 < p.NCols-1; c1++ {
		cc1 := p.NCols - c1 - 1
		for c2 := c1 + 1; c2 < p.NCols; c2++ {
			// Append one multiplier per truth-table row.
			cc2 := p.NCols - c2 - 1
			mults := make([]clp.Nonzero, 0, nr/2)
			for r := 0; r < nr; r++ {
				if (r>>cc1)&1 == 1 && (r>>cc2)&1 == 1 {
					mults = append(mults, clp.Nonzero{
						Index: r,
						Value: 1.0,
					})
				}
			}

			// Append the complete column to the matrix.
			mat.AppendColumn(mults)
		}
	}

	// Append an additional column of -1s for the constant term.  That is,
	// we want to constrain all valid rows to equal k, for some unknown k.
	// But because we can't have a variable on the right-hand side (as in
	// A+B+...+Z = k), we move the variable to the left-hand side (as in
	// A+B+...+Z-k = 0).
	mults := make([]clp.Nonzero, nr)
	for r := 0; r < nr; r++ {
		mults[r] = clp.Nonzero{
			Index: r,
			Value: -1.0,
		}
	}
	mat.AppendColumn(mults)

	// Append an additional column for the gap.  This is -1 for invalid
	// rows and 0 (omitted) for valid rows.
	mults = make([]clp.Nonzero, 0, nr)
	for r := 0; r < nr; r++ {
		if !isValid[r] {
			mults = append(mults, clp.Nonzero{
				Index: r,
				Value: -1.0,
			})
		}
	}
	mat.AppendColumn(mults)

	// Bound valid rows to [0, 0] and invalid rows to [epsilon, infinity].
	rb := make([]clp.Bounds, nr)
	positive := clp.Bounds{
		Lower: 1e-30, // Arbitrary small number, much larger than machine epsilon
		Upper: math.Inf(1),
	}
	for r, v := range isValid {
		if !v {
			rb[r] = positive
		}
	}

	// Bound each variable to its target range.
	cb := make([]clp.Bounds, nc+2)
	for c := range cb {
		switch {
		case c == nc+1:
			// Gap
			cb[c] = positive
		case c == nc:
			// Constant
			cb[c].Lower = math.Inf(-1)
			cb[c].Upper = math.Inf(1)
		case c < p.NCols:
			// Linear
			cb[c].Lower = p.MinL
			cb[c].Upper = p.MaxL
		default:
			// Quadratic
			cb[c].Lower = p.MinQ
			cb[c].Upper = p.MaxQ
		}
	}

	// The objective function is the gap, which we want to maximize.
	obj := make([]float64, nc+2)
	obj[nc+1] = 1.0

	// Solve the maximization problem.
	simp := clp.NewSimplex()
	simp.LoadProblem(mat, cb, obj, rb, nil)
	simp.SetOptimizationDirection(clp.Maximize)
	if simp.Primal(clp.NoValuesPass, clp.NoStartFinishOptions) != clp.Optimal {
		return false
	}
	soln := simp.PrimalColumnSolution()
	copy(q.Coeffs, soln)
	return true
}
