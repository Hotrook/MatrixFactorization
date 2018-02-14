# Sebastian Mroz 
# 221433 
function prettyPrintMatrix( A::Array{Float64, 2} )
    sizes = size(A)
    m = sizes[ 1 ]
    n = sizes[ 2 ]

    for i = 1 : m
        for j = 1 : n
            @printf("%7.1f ", A[ i, j ] )
        end
        println()
    end
    println()
end

function prettyPrintMatrix( A::CCSparseMatrix )
    prettyPrintMatrix( toDense( A ) )
end
