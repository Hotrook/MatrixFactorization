# Sebastian Mroz 
# 221433 
function columnCounts(A::CCSparseMatrix, parent::Array{Int64, 1}, post::Array{Int64, 1})
    m = A.m
    n = A.n
    AT = Sparse.transpose( A )

    delta = zeros(Int, n)
    colcount = zeros( Int, n )
    ancestor = collect(1:n)
    maxfirst = zeros( Int, n )
    prevleaf = zeros( Int, n )
    first = zeros( Int, n )

    initDeltaAndFirst!(post, parent, delta, first, n )

    for k = 1 : n
        j = post[ k ]
        if parent[ j ] != 0
            delta[ parent[ j ] ] -= 1
        end

        for p = AT.columnPointer[ j ] : (AT.columnPointer[ j+1 ] -1 )
            i = AT.rowIndex[ p ]
            q, jleaf = leaf!(i, j, first, maxfirst, prevleaf, ancestor )
            if jleaf >= 1
                delta[ j ] += 1
            end
            if jleaf == 2
                delta[ q ] -= 1
            end
        end

        if parent[ j ] != 0
            ancestor[ j ] = parent[ j ]
        end
    end

    for j = 1 : n
        if parent[ j ] != 0
            delta[ parent[ j ] ] += delta[ j ]
        end
    end

    return delta
end

function initDeltaAndFirst!(post::Array{Int64,1 }, parent::Array{Int64,1 },
    delta::Array{Int64,1 }, first::Array{Int64,1 }, n::Int64)
    for k = 1 : n
        j = post[ k ]
        if first[ j ] == 0
            delta[ j ] = 1
        end
        while j != 0 && first[ j ] == 0
            first[ j ] = k
            j = parent[ j ]
        end
    end
end

function leaf!(i, j::Int64, first, maxfirst, prevleaf, ancestor::Array{Int64,1})
    jleaf = 0

    if i <= j || first[ j ] <= maxfirst[ i ]
        return -1, -1
    end

    maxfirst[ i ] = first[ j ]
    jprev = prevleaf[ i ]
    prevleaf[ i ] = j
    if jprev == 0
        return i, 1
    end
    jleaf = 2

    q = jprev
    while q != ancestor[ q ]
        q = ancestor[ q ]
    end
    s = jprev
    while s != q
        sparent = ancestor[ s ]
        ancestor[ s ] = q
        s = sparent
    end
    return q, jleaf
end
