# Sebastian Mroz 
# 221433 
function fkeep( A::CCSparseMatrix, filter::Function )
    nz = 1

    for j = 1 : A.n
        p = A.columnPointer[ j  ]
        A.columnPointer[ j  ] = nz
        while p < A.columnPointer[ j + 1 ]
            if filter( A.rowIndex[ p ], j, A.value[ p ] )
                A.value[ nz ] = A.value[ p  ]
                A.rowIndex[ nz ] = A.rowIndex[ p ]
                nz += 1
            end
            p += 1
        end
    end
    A.columnPointer[ A.n+1 ] = nz
    A.nzmax = nz
    return nz
end
