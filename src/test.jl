include("Sparse.jl")
include("BenchmarkHelpers.jl")

using Sparse
using BenchmarkHelpers

in_directory = "testdata/matrices/positive_definite_matrices"
algorithm = "cholesky"
out_directory = "testresults/netlib/" * algorithm * "/"
probes = 10

initFunctions()
testMatricesFromDirectory( in_directory, out_directory, probes, algorithm )
