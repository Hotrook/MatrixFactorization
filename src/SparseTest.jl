# Sebastian Mroz 
# 221433 
using Sparse
using SuiteSparse
using Base.Test
import Base.transpose

function generateExampleMatrix()
    A = Array{Float64, 2}(3, 3)

    A[1,1] = 1.0
    A[2,1] = 0.0
    A[3,1] = 0.0
    A[1,2] = 2.0
    A[2,2] = 2.0
    A[3,2] = 2.0
    A[1,3] = 3.0
    A[2,3] = 0.0
    A[3,3] = 3.0

    return A
end

function equal( A::Sparse.CCSparseMatrix, B::Sparse.CCSparseMatrix )
    if A.nzmax == B.nzmax && A.n == B.n &&
        A.m == B.m &&
        A.columnPointer == B.columnPointer &&
        A.rowIndex == B.rowIndex &&
        A.value == B.value
        return true
    end
    return false
end

function makeSymmetric( A, n )
    for i = 1 : n, j = i : n
        A[ j, i ] = A[ i, j ]
    end
end

function readMatrixFromFile( filename::String )
    file = open( filename )
    lines = readlines( file )

    for j = 1 : 4
        lines[ j ] = lines[ j ][ 2 : ( length( lines[ j ] ) - 1) ]
        lines[ j ] = replace( lines[ j ] , ",", "")
    end

    columnPointer = map( x -> parse(Int64, x) , split(lines[ 1 ]) )
    rowIndex = map( x -> parse(Int64, x), split( lines[ 2 ] ) )
    value = map( x -> parse(Float64, x), split( lines[ 3 ] ) )
    expected = map( x -> parse(Int64, x) + 1, split( lines[ 4 ] ) )

    n = length(columnPointer ) - 1
    A = CCSparseMatrix( columnPointer[ n + 1 ] - 1, n, n )
    A.columnPointer= columnPointer
    A.rowIndex = rowIndex
    A.value = value

    close( file )

    return A, expected
end

@testset "SparseMatrix" begin
    maxSize = 10
    n = 20
    A = Sparse.CCSparseMatrix( 10, 20, 20)

    @test size(A.columnPointer)[1] == n+1
    @test size(A.rowIndex)[1] == maxSize
    @test size(A.rowIndex)[1] == maxSize
end

@testset "toSparse" begin
    A = generateExampleMatrix()
    B = Sparse.toSparse(A)

    @test B.columnPointer == [1,2,5,7]
    @test B.rowIndex == [1,1,2,3,1,3]
    @test B.value == [1.0,2.0,2.0,2.0,3.0,3.0]
end

@testset "cholesky" begin
    @testset "cholesky_without_pinv" begin
        expectedL = generateL(100, 0.3)
        A = Sparse.generatePosDefMatrix( expectedL )
        sparse_A = toSparse( A )

        sparse_L = cholesky_fact_without_analisys( sparse_A )
        L = toDense( sparse_L )

        @test expectedL == L
    end

    @testset "cholesky_with_pinv" begin

        for t = 1:9
            file = string("testdata/amdtestset", t, ".txt")
            css_a, e = readMatrixFromFile( file )
            expected = toDense( css_a )

            N, S = cholesky_fact( css_a )

            css_lt = Sparse.transpose( N.L )
            css_r = multiply( N.L, css_lt )
            css_actual = symperm( css_r, Sparse.pinv(S.pinv, 10) )

            actual = toDense( css_actual)
            makeSymmetric( actual, 10 )

            for i = 1:10, j = 1:10
                @test expected[ i, j ] ≈ actual[ i, j ]
            end
        end
    end

    @testset "cholesky_with_gen" begin
        n = 1000
        expected = generatePosDefMatrix(n, 0.01)
        css_a = toSparse(expected)

        N, S = cholesky_fact( css_a )

        css_lt = Sparse.transpose( N.L )
        css_r = multiply( N.L, css_lt )
        css_actual = symperm( css_r, Sparse.pinv(S.pinv, n) )

        actual = toDense( css_actual)
        makeSymmetric( actual, n )

        for i = 1:n, j = 1:n
            @test expected[ i, j ] ≈ actual[ i, j ]  atol=0.00000000001
        end
    end
end


@testset "help_functions" begin

    @testset "add" begin
        regularA = generateExampleMatrix()
        A = toSparse( regularA )
        B = toSparse( regularA )

        C = add(A, B, 1.0, 1.0)
        regularC = toDense( C )

        for i = 1:3, j = 1:3
            @test regularC[ i, j ] == regularA[ i, j] * 2
        end
    end

    @testset "multiply" begin
        A = zeros( Float64, 3, 3)
        B = zeros( Float64, 3, 3)
        expected = zeros( Float64,  3, 3)
        expected[ 1, 1 ] = 3.0

        for i = 1:3
            A[ 1, i ] = 1.0
            B[ i, 1 ] = 1.0
        end

        sparse_A = toSparse( A )
        sparse_B = toSparse( B )

        sparse_result = multiply(sparse_A, sparse_B )
        result = toDense( sparse_result )

        for i = 1:3, j = 1:3
            @test expected[ i, j ] == result[ i, j ]
        end
    end

    @testset "transpose" begin
        example = generateExampleMatrix()
        A = Sparse.toSparse( example )
        AT = Sparse.transpose( A )
        ATT = Sparse.transpose( AT )
        @test equal(A, ATT)

        M  = generatePosDefMatrix( 50, 0.5)
        A = Sparse.toSparse( M )
        AT = Sparse.transpose( A )
        MT = toDense( AT )

        @test Base.transpose( M ) == MT

    end

    @testset "pinv" begin
        @testset "first-case" begin
            P = [1,2,3,4,5]

            actual = Sparse.pinv( P, 5 )

            @test P == actual
        end

        @testset "second-case" begin
            P = [2, 3, 1, 5, 4]
            expected = [3, 1, 2, 5, 4]

            actual = Sparse.pinv( P, 5 )

            @test expected == actual
        end
    end
end

@testset "amd" begin

   @testset "5size" begin
        A = CCSparseMatrix( 11, 5, 5 )
        A.columnPointer = [1, 4, 5, 6, 9, 12]
        A.rowIndex = [1, 4, 5, 2,3,1, 4, 5,1, 4, 5,]
        A.value = [2704.0, 4940.0, 3380.0,8836.0,8100.0,4940.0, 13381.0, 6175.0, 3380.0, 6175.0, 6250.0]

        P = amd( 1, A )

        @test P == [ 2, 3, 1, 4, 5, 6 ]
    end

   @testset "amd-file" begin
        for i = 1 : 9
            filename = string("testdata/amdtestset" , string( i ), ".txt")
            A, expected = readMatrixFromFile( filename )

            P = amd( 1, A )

            @test P == expected
        end
    end

end

@testset "symperm" begin
    A = zeros(Float64, 3, 3)
    expected = zeros(Float64, 3, 3)
    expected[ 1, 1 ] = 4.0
    expected[ 1, 2 ] = 6.0
    expected[ 2, 2 ] = 9.0
    expected[ 1, 3 ] = 2.0
    expected[ 2, 3 ] = 3.0
    expected[ 3, 3 ] = 1.0
    P = [ 3, 1, 2 ]
    for i = 1 : 3, j = 1 : 3
        A[ i, j ] = i * j
    end
    sparse_A = toSparse( A )

    sparse_C = symperm(sparse_A, P )

    C = toDense( sparse_C )

    @test expected == C
end

@testset "lu" begin

    for i = 1 : 9
        file = string("testdata/amdtestset", i, ".txt")
        ccs_a, r = readMatrixFromFile(file)

        res, S = lu_fact( ccs_a )

        ccs_lu = multiply(res.L, res.U)
        ccs_expected = Sparse.permute( ccs_a, res.pinv, S.q)

        lu = toDense( ccs_lu )
        expected = toDense( ccs_expected )

        for i = 1 : 10, j = 1:10
            @test lu[ i, j ] ≈ expected[ i, j ] atol=0.00000000001
        end
    end

end
