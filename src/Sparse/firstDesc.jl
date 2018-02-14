# Sebastian Mroz 
# 221433 
function firstdesc(parent, post::Array{Int64, 1})
    n = size(parent)[ 1 ]
    first = zeros(Int64, n)
    level = fill(-1, n)

    for k = 1:n
        i = post[ k ]
        len = 0

        r = i
        while r != 0 && first[ r ] == 0
            first[ r ] = k
            r = parent[ r ]
            len = len + 1
        end

        if r == 0
            len = len - 1
        else
            len = len + level[ r ]
        end

        s = i
        while s != r
            level[ s ] = len
            len = len - 1
            s = parent[ s ]
        end
    end
    return first, level
end
