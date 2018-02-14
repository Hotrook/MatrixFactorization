# Sebastian Mroz 
# 221433 
function eliminationTree(A::CCSparseMatrix)
    parent = zeros(Int, A.n)
    ancestor = zeros(Int, A.n)

    for root = 1:A.n

        for p = (A.columnPointer[root]) : (A.columnPointer[root+1] - 1 )
            node = A.rowIndex[ p ]
            while node != 0 && node < root
                inext = ancestor[ node ]
                ancestor[ node ] = root
                if inext == 0
                    parent[ node ] = root
                end
                node = inext
            end

        end
    end

    return parent
end
