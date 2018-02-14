# Sebastian Mroz 
# 221433 
function ereach(A::CCSparseMatrix, k::Int, parent::Array{Int, 1}, stack::Array{Int, 1}, markspace::Array{Int, 1})
    top = A.n+1
    len = 0

    MARK(markspace, k)
    for p = A.columnPointer[k]:(A.columnPointer[k+1]-1)
        i = A.rowIndex[p]
        if i > k
            continue
        end

        len = 0
        while !MARKED(markspace, i)
            len += 1
            stack[ len ] = i
            MARK(markspace, i)
            i = parent[ i ]
        end
        while len >= 1
            top -= 1
            stack[ top ] = stack[ len ]
            len -= 1
        end
    end
    for p = top : A.n
        MARK( markspace, stack[ p ])
    end
    MARK( markspace, k )
    return top
end
