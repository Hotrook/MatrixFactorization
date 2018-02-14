# Sebastian Mroz 
# 221433 
function wclear( mark::Int64, lemax::Int64, w::Array{Int64, 1}, n::Int64)
    if mark < 2 || ( mark + lemax < 0 )
        for k = 1 : n
            if w[ k ] != 0
                w[ k ] = 1
                mark = 2
            end
        end
    end
    return mark
end
