# Sebastian Mroz 
# 221433 
function add( A::CCSparseMatrix, B::CCSparseMatrix, alpha::Float64, beta::Float64 )
    m = A.m
    anz = A.columnPointer[ A.n+1 ]-1
    n =  B.n
    bnz = B.columnPointer[ n+1 ]-1
    w = zeros( Int64, m+1 )

    C = CCSparseMatrix(anz + bnz, m, n )
    value = zeros( Float64, m+1 )

    nz = 1
    for j = 1:n
        C.columnPointer[ j ] = nz
        nz = scatter( A, j, alpha, w, value, j + 1, C, nz )
        nz = scatter( B, j, beta, w, value, j + 1, C, nz )
        for p = C.columnPointer[ j ]:( nz - 1 )
            C.value[ p ] = value[ C.rowIndex[ p ] ]
        end
    end

    C.columnPointer[ n+1 ] = nz
    return C
end


function addWithoutValues( A::CCSparseMatrix, B::CCSparseMatrix, alpha::Float64, beta::Float64 )
    m = A.m
    anz = A.columnPointer[ A.n+1 ]-1
    n =  B.n
    bnz = B.columnPointer[ n+1 ]-1
    w = zeros( Int64, m+1 )

    C = CCSparseMatrix(anz + bnz, m, n )
    value = zeros( Float64, m+1 )

    nz = 1
    for j = 1:n
        C.columnPointer[ j ] = nz
        nz = scatterWithoutValues( A, j, alpha, w, value, j + 1, C, nz )
        nz = scatterWithoutValues( B, j, beta, w, value, j + 1, C, nz )
    end

    C.columnPointer[ n+1 ] = nz
    return C
end
