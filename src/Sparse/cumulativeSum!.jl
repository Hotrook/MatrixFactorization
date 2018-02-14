# Sebastian Mroz 
# 221433 
function cumulativeSum!(p, c::Array{Int64, 1 }, n::Int64)
    nz = 1
    for i = 1 : n
        p[ i ] = nz
        nz = nz + c[ i ]
        c[ i ] = p[ i ]
    end
    p[ n+1 ] = nz
    return nz
end
