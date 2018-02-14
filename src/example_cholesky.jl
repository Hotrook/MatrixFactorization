# Sebastian Mroz
# 221433
using Sparse

n = 10
density = 0.4


A = generatePosDefMatrix(n, density)
sparse_a = toSparse( A )
    println("Matrix A:")
    prettyPrintMatrix( A )


# Przykład użycia funkcji cholesky_fact
N, S = cholesky_fact( sparse_a )
sparse_lt = transpose( N.L )
    println("Matrix L:")
    prettyPrintMatrix( N.L )
    println("Matrix LT:")
    prettyPrintMatrix( sparse_lt )

sparse_a2 = multiply( N.L, sparse_lt )
    println("Matrix L * LT:")
    prettyPrintMatrix( sparse_a2 )

sparse_pap = symperm(sparse_a, S.pinv)
    println("Matrix P * A * PT:")
    prettyPrintMatrix( sparse_pap )
