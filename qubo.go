// This file tries to find the QUBO parameters that implement a given
// truth table.

package main

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

// A QUBO represents a QUBO matrix whose values we're solving for.
type QUBO struct {
	Params *Parameters // Pointer to global program parameters
	Coeffs []float64   // List of linear followed by quadratic coefficients
}

// AllPossibleColumns returns a matrix containing all possible columns
// for a given number of rows that contain only 0s and 1s.
func AllPossibleColumns(nr int) *mat.Dense {
	nc := 1 << nr
	all := mat.NewDense(nr, nc, make([]float64, nr*nc))
	for r := 0; r < nr; r++ {
		for c := 0; c < nc; c++ {
			if c&(1<<r) != 0 {
				all.Set(r, c, 1.0)
			}
		}
	}
	return all
}

// Evaluate computes the badness of a set of coefficients.
func (q QUBO) Evaluate() (float64, error) {
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

	// Temporary
	{
		for r := 0; r < n; r++ {
			notify.Printf("Q[%d] = [%4.1f %4.1f %4.1f]",
				r, Q.At(r, 0), Q.At(r, 1), Q.At(r, 2))
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
	notify.Printf("Vals = %v", vals) // Temporary

	// Find the minimum output across all inputs.
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
