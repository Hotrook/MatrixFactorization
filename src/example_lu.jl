# Sebastian Mroz 
# 221433 
using Sparse

n = 10
density = 0.4


# Przykład użycia funkcji lu_fact
A = generatePosDefMatrix(n, density)
sparse_a = toSparse( A )
    println("Matrix A:")
    prettyPrintMatrix( sparse_a )

N, S = lu_fact( sparse_a )
sparse_lt = transpose( N.L )
    println("Matrix L:")
    prettyPrintMatrix( N.L )
    println("Matrix U:")
    prettyPrintMatrix( N.U )

sparse_a2 = multiply( N.L, N.U )
    println("Matrix L * U:")
    prettyPrintMatrix( sparse_a2 )

sparse_paq = permute(sparse_a, N.pinv, S.q)
    println("Matrix P * A * Q:")
    prettyPrintMatrix( sparse_paq )
