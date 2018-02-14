# Sebastian Mroz 
# 221433 
type Symbolic
    C::CCSparseMatrix
    pinv::Array{Int64, 1}
    q::Array{Int64, 1}
    parent::Array{Int64, 1}
    columnPointers::Array{Int64, 1}
    lnz::Float64
    unz::Float64

    Symbolic() = new()
end
