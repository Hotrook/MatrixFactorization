# Sebastian Mroz 
# 221433 
function isNotZero(i, j, aij)
    return aij != 0
end

function multiply(A::CCSparseMatrix, B::CCSparseMatrix)

    nz = 1
    m = A.m
    anz = A.columnPointer[ A.n+1 ]
    n = B.n
    bnz = B.columnPointer[ n+1 ]
    C = CCSparseMatrix( anz + bnz, m, n )

    w = zeros( Int64, m )
    x = zeros( Float64, m )

    for j = 1 : n
        if nz + m > C.nzmax
            resize!(C , C.nzmax + m  )
        end
        C.columnPointer[ j ] = nz
        for p = (B.columnPointer[ j ]) : (B.columnPointer[ j+1 ] - 1 )
            nz = scatter( A, B.rowIndex[ p ], B.value[ p ], w, x, j + 1, C, nz )
        end
        for p = C.columnPointer[ j ] : (nz-1)
            C.value[ p ] = x[ C.rowIndex[ p ] ]
        end
    end

    C.columnPointer[ n+1 ] = nz
    # fkeep( C, isNotZero )
    return C
end


function multiplyWithoutValues(A::CCSparseMatrix, B::CCSparseMatrix)

    nz = 1
    m = A.m
    anz = A.columnPointer[ A.n+1 ]
    n = B.n
    bnz = B.columnPointer[ n+1 ]
    C = CCSparseMatrix( anz + bnz, m, n, -1 )

    w = zeros( Int64, m )
    x = zeros( Float64, 0 )

    for j = 1 : n
        if nz + m > C.nzmax
            resize!(C.rowIndex , C.nzmax + m  )
            C.nzmax = C.nzmax + m
        end
        C.columnPointer[ j ] = nz
        for p = (B.columnPointer[ j ]) : (B.columnPointer[ j+1 ] - 1 )
            nz = scatterWithoutValues( A, B.rowIndex[ p ], 1.0, w, x, j + 1, C, nz )
        end
    end

    C.columnPointer[ n+1 ] = nz
    return C
end
