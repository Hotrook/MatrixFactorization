# Sebastian Mroz 
# 221433 
function symbolic_lu( order::Int64, A::CCSparseMatrix )
    S = Symbolic()
    n = A.n
    S.q = amd( order, A )
    S.unz = 4 * A.columnPointer[ n + 1 ] + n
    S.lnz = S.unz
    return S
end
