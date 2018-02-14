# Sebastian Mroz 
# 221433 
function spsolve( G::CCSparseMatrix, B::CCSparseMatrix, k::Int64, xi::Array{Int64, 1},
     x::Array{Float64, 1}, pinv::Array{Int64, 1}, lo::Int64, stack::Array{Int64, 1})

    n = G.n
    top = reach(G, B, k, xi, pinv, stack)
    p = q = 0

    for p = top : n
        x[ xi[ p ] ] = 0
    end
    for p = B.columnPointer[ k ] : (B.columnPointer[ k + 1 ] - 1 )
        x[ B.rowIndex[ p ] ] = B.value[ p ]
    end
    for px = top : n
        j = xi[ px ]
        J = pinv[ j ]
        if J < 0
            continue
        end
        index = G.columnPointer[ J ]
        if lo == 0
            index = G.columnPointer[ J + 1 ] - 1
        end

        x[ j ] /= G.value[ index ]
        if lo == 1
            p = G.columnPointer[ J ] + 1
            q = G.columnPointer[ J + 1 ]
        else
            p = G.columnPointer[ J ]
            q = G.columnPointer[ J + 1 ] - 1
        end
        for i = p : ( q - 1 )
            x[ G.rowIndex[ i ] ] -= G.value[ i ] * x[ j ]
        end
    end

    return top
end
