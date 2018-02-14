# Sebastian Mroz 
# 221433 
function generatePosDefMatrix(n::Int64)
    generatePosDefMatrix(n, 1.0)
end

function generatePosDefMatrix(n::Int64, coverage::Float64)
    A = rand( n, n)
    A += Base.transpose(A)
    A /= 2

    for i = 1:n, j = i:n
        if rand() > coverage
            A[ i, j ] = A[ j, i ] =0.0
        end
    end

    A += eye(n) * n

    A
end

function generatePosDefMatrix( L::Array{Float64, 2} )
    L * Base.transpose(L)
end

function generateL( size, coverage )
    L = zeros( size, size )
    for i = 1 : size
        for j = 1 : i
            if rand() < coverage  || i == j
                x = rand(1:10)
                L[ i, j ] = x
            end
        end
    end
    L
end

function generateRandomMatrix( size, coverage )
    A = rand( size, size)

    for i = 1:size, j = i:size
        if rand() > coverage && i != j
            A[ i, j ] = A[ j, i ] = 0.0
        end
    end
    A
end
