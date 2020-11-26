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

// NearZero defines a value close enough to zero to be considered a zero gap.
const NearZero = 1e-5

// A QUBO represents a QUBO matrix whose values we're solving for.
type QUBO struct {
	Params *Parameters // Pointer to global program parameters
	Coeffs []float64   // List of linear followed by quadratic coefficients
}

// QUBOFactory produces all QUBOs we care about in turn.  It returns in turn
// each QUBO with coefficients of the form n/c for n∈[-c, c] and c being the
// number of coefficients.
func QUBOFactory(p *Parameters) chan *QUBO {
	// Create a QUBO channel on which to send complete QUBOs.
	qch := make(chan *QUBO, 16)

	// Produce QUBOs in the background.
	go func() {
		// Iterate over increasing values of c on the hypothesis that
		// small values will suffice in most cases.
		nc := p.NCols               // Number of columns
		ncfs := (nc * (nc + 1)) / 2 // Number of coefficience
		icfs := make([]int, ncfs)
		for c := 1; c <= ncfs; c++ {
			// Create a channel on which to receive coefficient
			// lists.
			cch := make(chan []int, nc)

			// Create coefficient lists in the background.
			go func() {
				allCoeffs(icfs, 0, c, cch)
			}()

			// In the foreground (of a background goroutine), pack
			// coefficient lists into QUBOs for the caller's
			// convenience.
			for cf := range cch {
				// Divide each integer by nc to produce a
				// floating-point coefficient.
				cfs := make([]float64, ncfs)
				for i, c := range cf {
					cfs[i] = float64(c) / float64(ncfs)
				}

				// Create a QUBO and send it down the channel.
				qch <- &QUBO{
					Params: p,
					Coeffs: cfs,
				}
			}
		}
		close(qch)
	}()

	// Return a channel from which QUBOs can be read.
	return qch
}

// allCoeffs returns all possible slices with coefficients in the range
// [-c, c], starting with values of smallest magnitude.
func allCoeffs(cfs []int, i, c int, ch chan []int) {
	// Recursively fill in elements of cfs until none remain.
	if i == len(cfs) {
		// Ensure that at least one coefficient is equal to s.
		ok := false
		for _, cf := range cfs {
			if cf == c || cf == -c {
				ok = true
				break
			}
		}
		if !ok {
			return
		}

		// Send the completed coefficient list into the channel.
		coeffs := make([]int, i)
		copy(coeffs, cfs)
		ch <- coeffs
		return
	}
	for j := 0; j < 2*c+1; j++ {
		cfs[i] = j*(j%2) - j/2 // {0, 1, -1, 2, -2, ...}
		allCoeffs(cfs, i+1, c, ch)
	}
	if i == 0 {
		close(ch)
	}
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

// SelectValidRows selects at most one row in each batch of 2^a (with
// a ancillae) to treat as valid.
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

// ComputeGap returns the difference between the minimum invalid and maximum
// valid rows.  A positive difference indicates a correct solution; a
// non-positive difference indicates an incorrect solution.
func ComputeGap(vals []float64, isValid []bool) float64 {
	minInvalid, maxValid := math.MaxFloat64, -math.MaxFloat64
	for r, v := range vals {
		if isValid[r] {
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
	// Consider a large number of coefficients in turn.
	qch := QUBOFactory(p)
	for q := range qch {
		// Evaluate x'*Q*x for all binary x to produce a
		// vector of values.
		vals := q.EvaluateAllInputs()

		// Given the vector of values, select the smallest of each
		// batch of equally acceptable valid rows (i.e., same truth
		// table but different ancillae) to treat as valid and mark the
		// remaining rows as invalid.
		isValid := q.SelectValidRows(vals)

		// Ensure that all valid rows have a value less than that of
		// any invalid row.
		gap := ComputeGap(vals, isValid)
		if gap >= 0.0 {
			// The valid and invalid rows overlap.  Continue with
			// the next QUBO.
			continue
		}

		// At this point we have a correct QUBO but one that ignored
		// the bounds imposed on the coefficients.  We run an LP solver
		// both to enforce coefficient bounds and to maximize the gap
		// between valid and invalid rows.
		if !q.LPReoptimize(isValid) {
			notify.Fatal("The LP solver failed to optimize the QUBO coefficients")
		}
		vals = q.EvaluateAllInputs()    // Recompute the values now that we have new coefficients.
		gap = ComputeGap(vals, isValid) // Recompute the gap, too.
		if gap <= 0.0 {
			// False alarm.  The LP solver thinks it solved the
			// problem, but this was in fact a bogus solution
			// caused by numerical imprecision.
			continue
		}
		return q, gap, vals, isValid
	}
	return nil, 0.0, nil, nil
}
