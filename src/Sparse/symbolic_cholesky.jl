# Sebastian Mroz 
# 221433 
function symbolic_cholesky( order::Int64, A::CCSparseMatrix)
    S = Symbolic()
    n = A.n

    P = amd( order, A )
    S.pinv = pinv( P, n )
    S.C = symperm( A, S.pinv )
    S.parent = eliminationTree( S.C )
    post = postordering(S.parent)
    c = columnCounts(S.C, S.parent, post)

    S.columnPointers = zeros( Int64, n+1 )
    S.unz = S.lnz = cumulativeSum!(S.columnPointers, c, n)

    if S.lnz < 0
        println("Something went wrong. S.lnz < 0")
    end

    return S
end
