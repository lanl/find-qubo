###################################
# Build find-qubo                 #
# By Scott Pakin <pakin@lanl.gov> #
###################################

# Use the following to enable MPI support.
GOFLAGS = -tags mpi

# The following are Open MPI specific.
export CGO_CFLAGS = $(shell mpicc --showme:compile)
export CGO_LDFLAGS = $(shell mpicc --showme:link)

SOURCES = \
	main.go \
	mpi.go \
	no-mpi.go \
	params.go \
	qubo.go	\
	ttable.go

all: find-qubo

find-qubo: $(SOURCES)
	go build $(GOFLAGS) -o find-qubo

check: find-qubo
	go vet
	golint

clean:
	$(RM) find-qubo

.PHONY: all setenv check clean
