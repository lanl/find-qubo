// This file tries to find the QUBO parameters that implement a given
// truth table.

package main

import (
	"fmt"
	"math"
	"math/rand"
	"os"
	"strings"
	"time"

	"github.com/lanl/clp"
	"github.com/schollz/progressbar/v3"
	"gonum.org/v1/gonum/mat"
)

// A QUBO represents a QUBO matrix whose values we're solving for.
type QUBO struct {
	Params *Parameters // Pointer to global program parameters
	Coeffs []float64   // List of linear followed by quadratic coefficients
}

// NewQUBO initializes and return a new QUBO structure.
func NewQUBO(p *Parameters) *QUBO {
	nc := p.NCols
	return &QUBO{
		Params: p,
		Coeffs: make([]float64, (nc*(nc+1))/2),
	}
}

// ProduceTruthTable randomly generates a truth table from an existing one.
// All rows above a high-water mark are set to "invalid" (maximum flexibility
// in solving).  All rows below a low-water mark are left unmodified.  All rows
// in between are randomized.
// preserving a given number of rows's valid state and randomizing the rest.
func ProduceTruthTable(p *Parameters, rng *rand.Rand, oldTT TruthTable, lo, hi int) TruthTable {
	tt := oldTT.Copy() // Truth table to modify and return
	nr := 1 << p.NCols // Number of rows in the truth table
	br := 1 << p.NAnc  // Number of rows in a block
	for r := lo; r < hi; r += br {
		// Select at most one row per block to consider valid.
		if !p.TT[r] {
			// Entire block of rows is invalid: Leave it alone.
			continue
		}
		vOfs := rng.Intn(br) // Select a single, random row in the block to keep valid.
		for ofs := 0; ofs < br; ofs++ {
			tt[r+ofs] = ofs == vOfs
		}
	}
	for r := hi; r < nr; r++ {
		tt[r] = false
	}
	return tt
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

// LPSolve uses a linear-programming algorithm to re-optimize a QUBO's
// coefficients in search of perfect balance of all valid rows.  This method
// returns a success code.
func (q *QUBO) LPSolve(isValid TruthTable) bool {
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

// ComputeGap returns the difference between the minimum invalid and maximum
// valid rows.  A positive difference indicates a correct solution; a
// non-positive difference indicates an incorrect solution.
func ComputeGap(vals []float64, tt TruthTable) float64 {
	minInvalid, maxValid := math.MaxFloat64, -math.MaxFloat64
	for r, v := range vals {
		if tt[r] {
			maxValid = math.Max(maxValid, v)
		} else {
			minInvalid = math.Min(minInvalid, v)
		}
	}
	return minInvalid - maxValid
}

// OptimizeCoeffs tries to find the coefficients that best represent the given
// truth table.  It returns the QUBO, invalid-valid gap, row values, and row
// validity indicators.  The function aborts on error.
func OptimizeCoeffs(p *Parameters) (*QUBO, float64, []float64, []bool) {
	// Create a progress bar to show our progress.
	bar := progressbar.NewOptions(p.NRands,
		progressbar.OptionSetWriter(os.Stderr),
		progressbar.OptionFullWidth(),
		progressbar.OptionThrottle(time.Second),
		progressbar.OptionSetRenderBlankState(true),
		progressbar.OptionClearOnFinish())
	defer bar.Clear()

	// Work on increasingly large subsets of the truth table.
	rng := rand.New(rand.NewSource(time.Now().UTC().UnixNano()))
	nr := 1 << p.NCols // Number of rows in the truth table
	br := 1 << p.NAnc  // Number of rows in a block
	tries := 0         // Number of attempts [0, p.NRands-1]
	baseTT := p.TT     // Baseline truth table
	lo, hi := 0, br    // Start by working on a single block
	for tries < p.NRands {
		// Generate truth tables until we succeed.
		status.Printf("Working on [%d, %d] of %d rows", lo, hi, nr) // Temporary
		for tries < p.NRands {
			// Update the progress bar.
			tries++
			bar.Add(1)

			// Create a truth table.  Ensure it contains at least
			// one valid row.
			tt := ProduceTruthTable(p, rng, baseTT, lo, hi)
			allFalse := true
			for _, b := range tt {
				if b {
					allFalse = false
					break
				}
			}
			if allFalse {
				lo, hi = hi, hi+br
				continue
			}

			// Solve for coefficients representing the truth table.
			q := NewQUBO(p)
			if !q.LPSolve(tt) {
				continue
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
				continue
			}

			// We found a solution.  If it includes all rows, we're
			// done.  Otherwise, keep what we have and increase the
			// number of rows considered.
			if hi == nr {
				return q, gap, vals, tt
			}
			baseTT = tt
			lo, hi = hi, hi+br
			break
		}
	}
	return nil, 0.0, nil, nil // Return unsuccessfully.
}
