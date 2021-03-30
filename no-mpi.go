// This file provides MPI stubs for when MPI support is unavailable

// +build !mpi

package main

// MPIInit initializes MPI.
func MPIInit() {
}

// MPIFinalize finalizes MPI.
func MPIFinalize() {
}

// MPICommSize returns the number of ranks in MPI_COMM_WORLD.
func MPICommSize() int {
	return 1
}

// MPICommRank returns the caller's rank in MPI_COMM_WORLD.
func MPICommRank() int {
	return 0
}

// MPIBcastInt broadcasts a single integer from rank 0 to all other ranks.
func MPIBcastInt(n int) int {
	return n
}

// An MPIOp represents an MPI reduction operation.
type MPIOp int

 // MPIOpSum indicates a reduction that produces a sum.
var MPIOpSum = MPIOp(1)

 // MPIOpMax indicates a reduction that returns the maximum value.
var MPIOpMax = MPIOp(2)

// MPIReduceInts reduces one or more integers to rank 0 from all ranks.
func MPIReduceInts(op MPIOp, in []int) []int {
	return in
}

// MPIAllreduceInts reduces one or more integers to all ranks from all ranks.
func MPIAllreduceInts(op MPIOp, in []int) []int {
	return in
}

// MPISendInts sends one or more integers to another rank.
func MPISendInts(vs []int, to int) {
}

// MPIRecvInts receives one or more integers frmo another rank.
func MPIRecvInts(vs []int, from int) {
}
