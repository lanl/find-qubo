// This file tries to find the QUBO parameters that implement a given
// truth table.

package main

import (
	"fmt"
	"math"
	"strings"

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
var nVarsToMatrix map[int]*mat.Dense

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

// AsOctaveMatrix returns the coefficients as a string that can be pasted into
// GNU Octave or MATLAB.
func (q *QUBO) AsOctaveMatrix() string {
	n := q.NVars
	i := q.NVars
	oct := make([]string, n)
	for r := 0; r < n; r++ {
		row := make([]string, n)
		for c := 0; c < n; c++ {
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

// LPSolve uses a linear-programming algorithm to re-optimize a QUBO's
// coefficients in search of perfect balance of all valid rows.  This method
// returns a success code.
func (q *QUBO) LPSolve(p *Parameters, tt TruthTable) bool {
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
		Lower: 1e-30, // Arbitrary small number, much larger than machine epsilon
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
		// False alarm.  The LP solver thinks it solved
		// the problem, but this was in fact a bogus
		// solution likely caused by numerical
		// imprecision.
		return 0, nil
	}
	return gap, vals
}
