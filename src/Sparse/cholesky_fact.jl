# Sebastian Mroz 
# 221433 
function cholesky_fact(A::CCSparseMatrix)
    S = symbolic_cholesky( 1, A )
    cholesky_fact(A, S), S
end

function cholesky_fact(A::CCSparseMatrix, S::Symbolic)
    N = Numeric()
    n = A.n
    parent = S.parent
    cp = S.columnPointers
    C = S.C

    c = zeros(Int, 2*n )
    x = zeros(Float64, n )
    L = CCSparseMatrix(cp[n+1], n, n )

    stack = zeros(Int64, A.n)
    markspace = zeros(Int, A.n)

    try
        N.pinv = S.pinv
    catch
    end

    for k = 1 : n
        L.columnPointer[ k ] = c[ k ] = cp[ k ]
    end

    for k = 1 : n
        stack_start = ereach( C, k, parent, stack, markspace)
        x[ k ] = 0

        for p = C.columnPointer[ k ] : (C.columnPointer[ k + 1 ] -1 )
            if C.rowIndex[ p ] <= k
                x[ C.rowIndex[ p ] ] = C.value[ p ]
            end
        end

        d = x[ k ]
        x[ k ] = 0

        for top = stack_start : n
            i = stack[ top ]
            lki = x[ i ] / L.value[ L.columnPointer[ i ] ]
            x[ i ] = 0
            for p = (L.columnPointer[ i ] + 1) : (c[ i ] - 1 )
                x[ L.rowIndex[ p ] ] -= L.value[ p ] * lki
            end

            d -= lki * lki
            p = c[ i ]
            c[ i ] += 1
            L.rowIndex[ p ] = k
            L.value[ p ] = lki
        end

        if d <= 0
            println( "Non pos def matrix ")
            return
        end
        p = c[ k ]
        c[ k ] += 1
        L.rowIndex[ p ] = k
        L.value[ p ] = sqrt( d )
    end
    L.columnPointer[ n+1 ] = cp[ n+1 ]
    N.L = L

    N
end

function cholesky_fact_without_analisys(A::CCSparseMatrix)
    n = A.n
    S = Symbolic()

    S.parent = eliminationTree( A )
    post = postordering( S.parent )
    colCounts = columnCounts(A, S.parent, post)

    cp = zeros(Int64, n+1 )
    cumulativeSum!(cp, colCounts, n)
    S.columnPointers = cp
    S.C = A

    N = cholesky_fact(A, S)
    N.L
end
