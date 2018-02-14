# Sebastian Mroz 
# 221433 
function dfs(j::Int64, G::CCSparseMatrix, top::Int64, xi::Array{Int64, 1}, pinv::Array{Int64, 1}, pstack::Array{Int64, 1})
    xi[ 1 ] = j
    head = 1
    while head >= 1
        j = xi[ head ]
        jnew = pinv[ j ]

        if !MARKED( G.columnPointer, j )
            MARK( G.columnPointer, j )
            if jnew < 0
                pstack[ head ] = 0
            else
                pstack[ head ] = UNFLIP(G.columnPointer[ jnew ])
            end
        end

        done = true
        p2 = 0
        if jnew >= 1
            p2 = UNFLIP(G.columnPointer[ jnew + 1 ] )
        end

        for p = pstack[ head ] : (p2-1)
            i = G.rowIndex[ p ]
            if MARKED(G.columnPointer, i )
                continue
            end
            pstack[ head ] = p
            head += 1
            xi[ head ] = i
            done = false
            break
        end
        if done
            head -= 1
            top -= 1
            xi[ top ] = j
        end

    end
    return top
end
