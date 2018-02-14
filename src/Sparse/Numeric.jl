# Sebastian Mroz 
# 221433 
type Numeric
    L::CCSparseMatrix
    U::CCSparseMatrix
    pinv::Array{Int64, 1}

    Numeric() = new()
end
