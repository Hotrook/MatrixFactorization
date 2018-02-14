include("Sparse.jl")
include("BenchmarkHelpers.jl")

using BenchmarkHelpers
using Sparse
using SuiteSparse
using Base.Test
import Base.transpose

n = 2000
probes = 100

coverage_start = 0.001
coverage_step = 0.001
coverage_stop = 0.002

initFunctions()

prefix = "testresults/lu/test1test/"

enumeration_file = open(prefix * "enumeration.txt", "w")
symbolic_results_file = open(prefix * "symbolic_results.txt", "w")
ccsparse_results_file = open(prefix * "ccsparse_results.txt", "w")
julia_sparse_results_file = open(prefix * "julia_sparse_results.txt", "w")
julia_results_file = open(prefix * "julia_results.txt", "w")

for i = coverage_start : coverage_step : coverage_stop
    symbolic_results = Result()
    ccsparse_results = Result()
    julia_sparse_results = Result()
    julia_results = Result()

    for t = 1:probes
        @printf("%10.5lf %10d \r", i, t )
        A = generateRandomMatrix(n, i)
        cc_sparse_a = toSparse( A )
        julia_sparse_a = sparse( A )

        S, time, bytes, gctim, memallocs = @timed symbolic_lu(0, cc_sparse_a)
        addResults(symbolic_results, time, bytes, probes)

        _, time, bytes, gctim, memallocs = @timed lu_fact(cc_sparse_a, S, 1.0)
        addResults(ccsparse_results, time, bytes, probes)

        _, time, bytes, gctim, memallocs = @timed lufact(julia_sparse_a)
        addResults(julia_sparse_results, time, bytes, probes)

        _, time, bytes, gctim, memallocs = @timed lu(A)
        addResults(julia_results, time, bytes, probes)
    end

    toWrite = @sprintf("%10.5lf\n", i)
    write(enumeration_file, toWrite )

    writeResults( symbolic_results_file, symbolic_results)
    writeResults( ccsparse_results_file, ccsparse_results)
    writeResults( julia_sparse_results_file, julia_sparse_results)
    writeResults( julia_results_file, julia_results)

end

close(enumeration_file)
close(symbolic_results_file)
close(ccsparse_results_file)
close(julia_sparse_results_file)
close(julia_results_file)
