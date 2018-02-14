# Sebastian Mroz 
# 221433 
function resize!( A :: CCSparseMatrix, nzmax::Int64)
    resize!( A.rowIndex, nzmax )
    resize!( A.value, nzmax )
    A.nzmax = nzmax
end
