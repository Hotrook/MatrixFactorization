# Sebastian Mroz 
# 221433 
# Sebastian Mroz 
# 221433 
type CCSparseMatrix
    nzmax::Int
    m::Int
    n::Int
    columnPointer::Array{Int}
    rowIndex::Array{Int}
    value::Array{Float64}

    function CCSparseMatrix(nzmax::Int, m::Int, n::Int)
        new(nzmax, m, n, zeros(Int, n+1), zeros(Int, nzmax), zeros( Float64,nzmax))
    end

    function CCSparseMatrix(nzmax::Int, m::Int, n::Int, value::Int)
        new(nzmax, m, n, zeros(Int, n+1), zeros(Int, nzmax),Array{Float64,1}())
    end
end
