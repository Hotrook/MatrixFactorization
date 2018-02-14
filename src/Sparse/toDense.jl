# Sebastian Mroz 
# 221433 
function toDense( A::CCSparseMatrix )
    M = zeros( A.m, A.n )

    for j = 1:A.n
        for v = A.columnPointer[ j ] : ( A.columnPointer[ j + 1 ] - 1 )
            i = A.rowIndex[ v ]
            M[ i, j ] = A.value[ v ]
        end
    end

    return M
end
