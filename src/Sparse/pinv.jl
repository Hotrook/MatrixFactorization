# Sebastian Mroz 
# 221433 
function pinv( permutation::Array{Int64, 1}, n::Int64)
    pinv = zeros(Int64, n)
    for k = 1 : n
        pinv[ permutation[ k ] ] = k
    end
    return pinv
end
