# Sebastian Mroz 
# 221433 
function transpose(A::CCSparseMatrix)
    AT = CCSparseMatrix( A.nzmax, A.n, A.m )
    workspace = zeros( Int, A.m )

    for p = 1 : (A.columnPointer[ A.n+1 ] - 1)
        workspace[ A.rowIndex[ p ] ] = workspace[A.rowIndex[p]] + 1
    end
    cumulativeSum!(AT.columnPointer, workspace, A.m )

    for j = 1 : A.n
        for p = A.columnPointer[ j ] : (A.columnPointer[ j+1 ] - 1)
            q = workspace[ A.rowIndex[ p ] ]
            workspace[ A.rowIndex[ p ] ]  = workspace[ A.rowIndex[ p ] ] + 1
            AT.rowIndex[ q ] = j
            AT.value[ q ] = A.value[ p ]
        end
    end

    return AT
end

function transposeWithoutValues(A::CCSparseMatrix)
    AT = CCSparseMatrix( A.nzmax, A.n, A.m, -1 )
    workspace = zeros( Int, A.m )

    for p = 1 : (A.columnPointer[ A.n+1 ] - 1)
        workspace[ A.rowIndex[ p ] ] = workspace[A.rowIndex[p]] + 1
    end
    cumulativeSum!(AT.columnPointer, workspace, A.m )

    for j = 1 : A.n
        for p = A.columnPointer[ j ] : (A.columnPointer[ j+1 ] - 1)
            q = workspace[ A.rowIndex[ p ] ]
            workspace[ A.rowIndex[ p ] ]  = workspace[ A.rowIndex[ p ] ] + 1
            AT.rowIndex[ q ] = j
        end
    end

    return AT
end
