# Sebastian Mroz 
# 221433 
function symperm( A::CCSparseMatrix, pinv::Array{Int64, 1} )
    n = A.n
    w = zeros( Int64, n )
    C = CCSparseMatrix(A.columnPointer[ n + 1 ], n, n )
    for j =  1 : n
        j2 = pinv[ j ]
        for p = A.columnPointer[ j ] : ( A.columnPointer[ j + 1 ] - 1)
            i = A.rowIndex[ p ]
            if i > j
                continue
            end
            i2 = pinv[ i ]
            w[ max( i2, j2 ) ] += 1
        end
    end
    cumulativeSum!(C.columnPointer, w, n )
    for j = 1 : n
        j2 = pinv[ j ]
        for p = A.columnPointer[ j ] : (A.columnPointer[ j +  1] - 1 )
            i = A.rowIndex[ p ]
            if i > j
                continue
            end
            i2 = pinv[ i ]
            q = w[ max( i2, j2 )]
            C.rowIndex[ q ] = min(i2, j2)
            w[  max( i2, j2 ) ] += 1
            C.value[ q ] = A.value[ p ]
        end
    end

    return C
end
