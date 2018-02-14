# Sebastian Mroz 
# 221433 
function reach(G::CCSparseMatrix, B::CCSparseMatrix, k::Int64, xi::Array{Int64, 1}, pinv::Array{Int64, 1}, stack::Array{Int64, 1})
    n = G.n
    top = n + 1
    for p = B.columnPointer[ k ] : (B.columnPointer[ k + 1 ] - 1 )
        if !MARKED(G.columnPointer, B.rowIndex[ p ] )
            top = dfs(B.rowIndex[ p ], G, top, xi, pinv, stack)
        end
    end
    for p = top : n
        MARK( G.columnPointer, xi[ p ] )
    end
    return top
end
