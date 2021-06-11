// This file provides minimal MPI support.  It assumes the MPI_ERRORS_ARE_FATAL
// error handler so no error codes are returned.

// +build mpi

package main

/*
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
*/
import "C"
import "unsafe"

// MPIInit initializes MPI.
func MPIInit() {
	var argc C.int
	C.MPI_Init(&argc, nil)
}

// MPIFinalize finalizes MPI.
func MPIFinalize() {
	C.MPI_Finalize()
}

// MPICommSize returns the number of ranks in MPI_COMM_WORLD.
func MPICommSize() int {
	var sz C.int
	C.MPI_Comm_size(C.MPI_COMM_WORLD, &sz)
	return int(sz)
}

// MPICommRank returns the caller's rank in MPI_COMM_WORLD.
func MPICommRank() int {
	var r C.int
	C.MPI_Comm_rank(C.MPI_COMM_WORLD, &r)
	return int(r)
}

// MPIBcastInt broadcasts a single integer from rank 0 to all other ranks.
func MPIBcastInt(n int) int {
	buf := C.long(n)
	C.MPI_Bcast(unsafe.Pointer(&buf), 1, C.MPI_LONG, 0, C.MPI_COMM_WORLD)
	return int(buf)
}

// An MPIOp represents an MPI reduction operation.
type MPIOp C.MPI_Op

 // MPIOpSum indicates a reduction that produces a sum.
var MPIOpSum = MPIOp(C.MPI_SUM)

 // MPIOpMax indicates a reduction that returns the maximum value.
var MPIOpMax = MPIOp(C.MPI_MAX)

// MPIReduceInts reduces one or more integers to rank 0 from all ranks.
func MPIReduceInts(op MPIOp, in []int) []int {
	n := len(in)
	sBuf := make([]C.long, n)
	rBuf := make([]C.long, n)
	for i, v := range in {
		sBuf[i] = C.long(v)
	}
	C.MPI_Reduce(unsafe.Pointer(&sBuf[0]), unsafe.Pointer(&rBuf[0]), C.int(n), C.MPI_LONG, C.MPI_Op(op), 0, C.MPI_COMM_WORLD)
	out := make([]int, n)
	for i, v := range rBuf {
		out[i] = int(v)
	}
	return out
}

// MPIAllreduceInts reduces one or more integers to all ranks from all ranks.
func MPIAllreduceInts(op MPIOp, in []int) []int {
	n := len(in)
	sBuf := make([]C.long, n)
	rBuf := make([]C.long, n)
	for i, v := range in {
		sBuf[i] = C.long(v)
	}
	C.MPI_Allreduce(unsafe.Pointer(&sBuf[0]), unsafe.Pointer(&rBuf[0]), C.int(n), C.MPI_LONG, C.MPI_Op(op), C.MPI_COMM_WORLD)
	out := make([]int, n)
	for i, v := range rBuf {
		out[i] = int(v)
	}
	return out
}

// MPISendInts sends one or more integers to another rank.
func MPISendInts(vs []int, to int) {
	n := len(vs)
	sBuf := make([]C.long, n)
	for i, v := range vs {
		sBuf[i] = C.long(v)
	}
	C.MPI_Send(unsafe.Pointer(&sBuf[0]), C.int(n), C.MPI_LONG, C.int(to), 0, C.MPI_COMM_WORLD)
}

// MPIRecvInts receives one or more integers frmo another rank.
func MPIRecvInts(vs []int, from int) {
	n := len(vs)
	rBuf := make([]C.long, n)
	for i, v := range vs {
		rBuf[i] = C.long(v)
	}
	var status C.MPI_Status
	C.MPI_Recv(unsafe.Pointer(&rBuf[0]), C.int(n), C.MPI_LONG, C.int(from), 0, C.MPI_COMM_WORLD, &status)
	for i, v := range rBuf {
		vs[i] = int(v)
	}
}
