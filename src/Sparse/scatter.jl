# Sebastian Mroz 
# 221433 
function scatter( A::CCSparseMatrix, j::Integer, beta::Float64,
        w::Array{Int64,1},
        x::Array{Float64, 1},
        mark::Integer,
        C::CCSparseMatrix,
        nz::Integer )
    for p = A.columnPointer[ j ] : ( A.columnPointer[ j+1 ] - 1 )
        i = A.rowIndex[ p ]
        if w[ i ] < mark
            w[ i ] = mark
            C.rowIndex[ nz ] = i
            nz += 1
            x[ i ] = beta * A.value[ p ]
        else
            x[ i ] += beta*( A.value[ p ])
        end
    end

    nz
end


function scatterWithoutValues( A::CCSparseMatrix, j::Integer, beta::Float64,
        w::Array{Int64,1},
        x::Array{Float64, 1},
        mark::Integer,
        C::CCSparseMatrix,
        nz::Integer )
    for p = A.columnPointer[ j ] : ( A.columnPointer[ j+1 ] - 1 )
        i = A.rowIndex[ p ]
        if w[ i ] < mark
            w[ i ] = mark
            C.rowIndex[ nz ] = i
            nz += 1
        end
    end

    nz
end
