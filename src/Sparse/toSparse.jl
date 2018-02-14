# Sebastian Mroz 
# 221433 
function toSparse(A::Array{Float64, 2})
    sizes = size(A)
    nonzeros = countnz(A)
    result = CCSparseMatrix(nonzeros, sizes[1], sizes[2])

    counter = 1
    result.columnPointer[1] = 1

    for j = 1:sizes[2]
        result.columnPointer[ j+1 ] = result.columnPointer[j]
        for i = 1:sizes[1]
            if A[i, j] != 0.0
                result.value[ counter ] = A[ i, j ]
                result.rowIndex[ counter ] = i
                result.columnPointer[ j+1 ] = result.columnPointer[ j+1 ] + 1
                counter = counter + 1
            end
        end
    end

    return result
end
