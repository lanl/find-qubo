// This file tries to find the QUBO parameters that implement a given
// truth table.

package main

import (
	"math"
	"math/bits"
	"sort"
	"sync/atomic"
	"time"

	"github.com/lanl/clp"
	"gonum.org/v1/gonum/mat"
)

// A QUBO represents a QUBO matrix whose values we're solving for.
type QUBO struct {
	NVars  int       // Number of variables
	Coeffs []float64 // List of linear followed by quadratic coefficients
}

// NewQUBO initializes and return a new QUBO structure that represents a given
// number of variables.
func NewQUBO(nv int) *QUBO {
	return &QUBO{
		NVars:  nv,
		Coeffs: make([]float64, (nv*(nv+1))/2),
	}
}

// nVarsToMatrix caches mappings from a variable count, nv, to a matrix of nv
// rows and 2^nv columns that represents all possible nv-bit vectors.
var nVarsToMatrix = make(map[int]*mat.Dense, 10)

// allPossibleInputs returns a binary matrix containing all possible columns
// for a given number of rows.
func (q *QUBO) allPossibleInputs(nv int) *mat.Dense {
	// Check if we've already computed an appropriate matrix.
	if m, found := nVarsToMatrix[nv]; found {
		return m
	}

	// Create and cache a new matrix.
	pnv := 1 << nv
	m := mat.NewDense(nv, pnv, make([]float64, nv*pnv))
	for r := 0; r < nv; r++ {
		for c := 0; c < pnv; c++ {
			if c&(1<<r) != 0 {
				m.Set(nv-r-1, c, 1.0)
			}
		}
	}
	nVarsToMatrix[nv] = m
	return m
}

// EvaluateAllInputs multiplies the QUBO by each input column in turn (i.e.,
// x'*Q*x for all x).
func (q *QUBO) EvaluateAllInputs() []float64 {
	// Convert the coefficients to an upper-triangular matrix.
	n := q.NVars // Number of columns and rows
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

	// Evaluate the matrix on all 2^n possible inputs for n variables.
	all := q.allPossibleInputs(n)
	nr := 1 << n
	vals := make([]float64, nr)
	for r := 0; r < nr; r++ {
		col := all.ColView(r)
		m := mat.NewDense(1, 1, []float64{0.0})
		m.Product(col.T(), Q, col)
		vals[r] = m.At(0, 0)
	}
	return vals
}

// LPSolve uses a linear-programming algorithm to re-optimize a QUBO's
// coefficients in search of perfect balance of all valid rows.  This method
// returns a success code.
func (q *QUBO) LPSolve(p *Parameters, tt TruthTable) bool {
	// Because CLP supports only "a >= b", not "a > b", we define
	// epsilon as a small number and solve "a >= b + epsilon".
	simp := clp.NewSimplex()
	simp.SetPrimalTolerance(p.Tolerance)
	epsilon := simp.PrimalTolerance()

	// Define multipliers for each linear term's coefficients.  Each
	// multiplier will be either 0.0 (omitted) or 1.0.
	mat := clp.NewPackedMatrix()
	nc := len(q.Coeffs) // Number of coefficients
	nr := tt.NRows      // Number of rows
	for c := 0; c < tt.NCols; c++ {
		// Append one multiplier per truth-table row.
		cc := q.NVars - c - 1
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
	for c1 := 0; c1 < q.NVars-1; c1++ {
		cc1 := q.NVars - c1 - 1
		for c2 := c1 + 1; c2 < q.NVars; c2++ {
			// Append one multiplier per truth-table row.
			cc2 := q.NVars - c2 - 1
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
		if !tt.TT[r] {
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
		Lower: epsilon,
		Upper: math.Inf(1),
	}
	for r, v := range tt.TT {
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
		case c < q.NVars:
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
	simp.LoadProblem(mat, cb, obj, rb, nil)
	simp.SetOptimizationDirection(clp.Maximize)
	if simp.Primal(clp.NoValuesPass, clp.NoStartFinishOptions) != clp.Optimal {
		return false
	}
	soln := simp.PrimalColumnSolution()
	copy(q.Coeffs, soln)
	return true
}

// ComputeGap returns the difference between the minimum invalid and maximum
// valid rows.  A positive difference indicates a correct solution; a
// non-positive difference indicates an incorrect solution.
func ComputeGap(vals []float64, tt TruthTable) float64 {
	minInvalid, maxValid := math.MaxFloat64, -math.MaxFloat64
	for r, v := range vals {
		if tt.TT[r] {
			maxValid = math.Max(maxValid, v)
		} else {
			minInvalid = math.Min(minInvalid, v)
		}
	}
	return minInvalid - maxValid
}

// trySolve attempts to solve for a QUBO's coefficients.  It returns the
// invalid-valid gap and row values.
func (q *QUBO) trySolve(p *Parameters, tt TruthTable) (float64, []float64) {
	// Run the solver.  Keep track of the number of invocations.
	atomic.AddUint64(&p.NumLPSolves, 1)
	if !q.LPSolve(p, tt) {
		return 0, nil
	}

	// Round the coefficients if asked to.
	rt := p.RoundTo
	if rt > 0.0 {
		for i, cf := range q.Coeffs {
			q.Coeffs[i] = math.Round(cf/rt) * rt
		}
	}

	// Compute the row values and gap.
	vals := q.EvaluateAllInputs()
	gap := ComputeGap(vals, tt)
	if gap <= 0.0 {
		// False alarm.  The LP solver thinks it solved the problem,
		// but this was in fact a bogus solution likely caused by
		// numerical imprecision.
		return 0, nil
	}
	return gap, vals
}

// allPossibleAncillae returns all length-N sequences of the numbers [0, k-1].
func allPossibleAncillae(n, k int) chan []int {
	ch := make(chan []int, 16)
	var tryNext func(work []int, idx int)
	tryNext = func(work []int, idx int) {
		// When the work slice is fully populated, send the final
		// sequence down the channel.
		if idx < 0 {
			final := make([]int, n)
			copy(final, work)
			ch <- final
			return
		}

		// Otherwise, try all possibilities for index idx, and
		// recursively try indices [0, idx-1].
		for v := 0; v < k; v++ {
			work[idx] = v
			tryNext(work, idx-1)
		}
	}
	go func() {
		tryNext(make([]int, n), n-1)
		close(ch)
	}()
	return ch
}

// bruteForceFindCoeffsWithAncillae uses a brute-force technique to solve for
// the QUBO coefficients given a specific number of ancillary variables to
// append to each row of the truth table.  The function returns the augmented
// truth table, the gap, the truth-table row values, the QUBO coefficients, and
// a success code.
func bruteForceFindCoeffsWithAncillae(p *Parameters, tt TruthTable, na int) (TruthTable, float64, []float64, []float64, bool) {
	// Create a truth table extended with ancillary variables and an
	// associated QUBO.
	ett := NewTruthTable(tt.NCols + na)
	q := NewQUBO(ett.NCols)

	// Lazily produce a list of all valid ancilla assignments to valid
	// truth-table rows.
	rows := tt.ValidRows()
	sortByOneBits(rows)
	naRows := 1 << na
	ch := allPossibleAncillae(len(rows), naRows)

	// Try in turn each possible set of ancillary variables.
	for ancs := range ch {
		// Set all ancillae according to ancs.
		ett.Clear()
		for i, r := range rows {
			rOfs := ancs[i]
			ett.TT[r*naRows+rOfs] = true
		}

		// See if the current configuration is solvable.
		gap, vals := q.trySolve(p, ett)
		if vals != nil {
			// Success!
			return ett, gap, vals, q.Coeffs, true
		}
	}
	return TruthTable{}, 0, nil, nil, false // Failure
}

// sortByOneBits sorts a given list of numbers by increasing number of 1 bits.
func sortByOneBits(nums []int) {
	// Precompute all 1-bit counts.
	n := len(nums)
	nBits := make(map[int]int, n)
	for _, v := range nums {
		nBits[v] = bits.OnesCount(uint(v))
	}

	// Sort the numbers.
	sort.SliceStable(nums, func(i, j int) bool {
		return nBits[nums[i]] < nBits[nums[j]]
	})
}

// findCoeffsWithAncillae solves for the QUBO coefficients given a specific
// number of ancillary variables to append to each row of the truth table.  The
// function returns the augmented truth table, the gap, the truth-table row
// values, the QUBO coefficients, and a success code.
func findCoeffsWithAncillae(p *Parameters, tt TruthTable, na int) (TruthTable, float64, []float64, []float64, bool) {
	// Create a truth table extended with ancillary variables and an
	// associated QUBO.
	ett := NewTruthTable(tt.NCols + na)
	q := NewQUBO(ett.NCols)

	// Sort the list of valid rows by increasing number of 1 bits.
	// Experimentation indicates that this may be a good approach for
	// finding a solution.
	rows := tt.ValidRows()
	sortByOneBits(rows)

	// Sort the numbers [0, naRows) by decreasing number of 1 bits.
	// Experimentation indicates that this may be a good approach for
	// finding a solution.
	naRows := 1 << na
	rowOfs := make([]int, naRows)
	for i := range rowOfs {
		rowOfs[i] = i
	}
	sortByOneBits(rowOfs)
	for i := 0; i < naRows/2; i++ {
		rowOfs[i], rowOfs[naRows-i-1] = rowOfs[naRows-i-1], rowOfs[i]
	}

	// Append one row at a time from the original truth table to the
	// extended truth table.
	var gap float64
	var vals []float64
RowLoop:
	for _, r := range rows {
		// Try in turn each possible set of ancillary variables.
		for _, rOfs := range rowOfs {
			// Try each bit pattern until one is solvable.
			ett.TT[r*naRows+rOfs] = true
			gap, vals = q.trySolve(p, ett)
			if vals == nil {
				// Failure: Undo the current ancilla pattern.
				ett.TT[r*naRows+rOfs] = false
			} else {
				// Success: Move on to the next row.
				continue RowLoop
			}
		}
		return TruthTable{}, 0, nil, nil, false // Failed to add the current row.
	}
	return ett, gap, vals, q.Coeffs, true // Success!
}

// FindCoefficients solves for the QUBO coefficients, adding ancillary
// variables as necessary.  It uses a greedy algorithm and does not guarantee
// that the number of ancillae is minimized.  The function returns the
// augmented truth table, the gap, the truth-table row values, the QUBO
// coefficients, and a success code.
func FindCoefficients(p *Parameters, tt TruthTable) (TruthTable, float64, []float64, []float64, bool) {
	// First check if the truth table is solvable without introducing any
	// ancillae.
	q := NewQUBO(tt.NCols)
	sTime := time.Now()
	gap, vals := q.trySolve(p, tt)
	eTime := time.Since(sTime)
	info.Printf("Performed %d LP solve in %d ms",
		p.NumLPSolves, eTime.Milliseconds())
	if vals != nil {
		return tt.Copy(), gap, vals, q.Coeffs, true
	}

	// Repeatedly increase the number of ancillae until we find a solution.
	for na := 1; na <= int(p.MaxAncillae); na++ {
		info.Printf("Increasing the number of ancillae to %d (%d total truth-table rows)", na, 1<<(tt.NCols+na))
		prevLP := p.NumLPSolves
		sTime := time.Now()
		var ett TruthTable   // Extended truth table
		var gap float64      // Minimal gap between valid and invalid rows
		var vals []float64   // QUBO evaluation on every possible input
		var coeffs []float64 // QUBO coefficients
		var ok bool          // Success code
		if p.BruteForce {
			ett, gap, vals, coeffs, ok = bruteForceFindCoeffsWithAncillae(p, tt, na)
		} else {
			ett, gap, vals, coeffs, ok = findCoeffsWithAncillae(p, tt, na)
		}
		eTime := time.Since(sTime)
		info.Printf("Performed %d LP solves in %d ms",
			p.NumLPSolves-prevLP, eTime.Milliseconds())
		if ok {
			return ett, gap, vals, coeffs, ok // Success!
		}
	}
	return TruthTable{}, 0.0, nil, nil, false // Failure
}
