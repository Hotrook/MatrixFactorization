# Sebastian Mroz 
# 221433 
function lu_fact( A::CCSparseMatrix )
    S = symbolic_lu( 2, A )
    lu_fact(A, S), S
end

function lu_fact( A::CCSparseMatrix, S::Symbolic )
    lu_fact(A, S, 1.0)
end

function lu_fact( A::CCSparseMatrix, S::Symbolic, tol::Float64)
    n = A.n
    unz = S.unz
    lnz = S.lnz
    N = Numeric()
    N.L = CCSparseMatrix( Int(lnz), n, n )
    N.U = CCSparseMatrix( Int(unz), n, n )
    N.pinv = pinv = fill(-1, n)
    x = zeros(Float64, n)
    a = t = pivot = 0.0
    col = 0
    lnz = unz = 1
    ipiv = 0
    xi = zeros(Int64, n)
    stack = zeros(Int64, n)
    i = top = p = 0

    for k = 1 : n
        N.L.columnPointer[ k ] = lnz
        N.U.columnPointer[ k ] = unz

        if lnz + n > N.L.nzmax
            resize!(N.L, lnz + n )
        end
        if unz + n > N.U.nzmax
            resize!(N.U, unz + n )
        end
        col = S.q[ k ]
        top = spsolve( N.L, A, col, xi, x, pinv, 1, stack )
        ipiv = -1
        a = -1.0
        for p = top : n
            i = xi[ p ]
            if pinv[ i ] < 0
                t = abs( x[ i ] )
                if t > a
                    a = t
                    ipiv = i
                end
            else
                N.U.rowIndex[ unz ] = pinv[ i ]
                N.U.value[ unz ] = x[ i ]
                unz += 1
            end
        end

        # println(a)
        if ipiv == -1 || a <= 0
            println("Something went wrong. This is a: $a\n")
            return
        end
        if pinv[ col ] < 0 && abs( x[ col ] ) >= a*tol
            ipiv = col
        end

        pivot = x[ ipiv ]
        N.U.rowIndex[ unz ] = k
        N.U.value[ unz ] = pivot
        unz += 1
        pinv[ ipiv ] = k
        N.L.rowIndex[ lnz ] = ipiv
        N.L.value[ lnz ] = 1
        lnz += 1

        for p = top : n
            i = xi[ p ]
            if pinv[ i ] < 0
                N.L.rowIndex[ lnz ] = i
                N.L.value[ lnz ] = x[ i ] / pivot
                lnz += 1
            end
            x[ i ] = 0
        end

    end

    N.L.columnPointer[ n + 1 ] = lnz
    N.U.columnPointer[ n + 1 ] = unz
    for p = 1 : (lnz - 1)
        N.L.rowIndex[ p ] = pinv[ N.L.rowIndex[ p ] ]
    end
    return N
end
