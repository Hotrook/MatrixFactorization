# Sebastian Mroz 
# 221433 
function permute(A::CCSparseMatrix, pinv::Array{Int64, 1}, q::Array{Int64, 1})
    n = A.n
    C = CCSparseMatrix(A.columnPointer[n+1], A.m, A.n)
    nz = 1

    for k = 1 : n
        C.columnPointer[ k ] = nz
        j = q[ k ]
        for t = A.columnPointer[ j ] : ( A.columnPointer[ j + 1 ] - 1 )
            C.value[ nz ] = A.value[ t ]
            C.rowIndex[ nz ] = pinv[ A.rowIndex[ t ] ]
            nz += 1
        end
        C.columnPointer[ n + 1 ] = nz
    end
    return C
end

function permute(pinv::Array{Int64, 1}, A::CCSparseMatrix)
    q = collect(1 :( length(pinv) ) )
    return permute( A, pinv, q )
end
