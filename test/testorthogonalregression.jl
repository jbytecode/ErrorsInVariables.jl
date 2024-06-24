using Test
using ErrorsInVariables
import Random

@testset "Orthogonal Regression" verbose = true begin

    @testset "R data example" begin

        eps = 0.0001

        X = Float64[1 1.0; 1 0.6; 1 1.2; 1 1.4; 1 0.2]
        y = Float64[0.5, 0.3, 0.7, 1.0, 0.2]

        result = orthogonal_regression(X, y)

        @test result.converged

        @test isapprox(result.betas[1], -0.03328271, atol=eps)
        @test isapprox(result.betas[2], 0.65145762, atol=eps)
    end



    @testset "Simple Model" begin

        rng = Random.MersenneTwister(1234)

        n = 30

        deltax = randn(rng, n) * sqrt(3.0)

        cleanx = randn(rng, n) * sqrt(7.0)

        e = randn(rng, n) * sqrt(5.0)

        y = 20.0 .+ 10.0 .* cleanx .+ e

        dirtyx = cleanx .+ deltax

        Xd = hcat(ones(n), dirtyx)
        Xc = hcat(ones(n), cleanx)
        dirtybetas = Xd \ y
        cleanbetas = Xc \ y

        result = orthogonal_regression(Xd, y, initialbetas=dirtybetas)

        dist_clean_and_real = euclidean(cleanbetas, [20.0, 10.0])
        dist_dirty_and_real = euclidean(dirtybetas, [20.0, 10.0])
        dist_eive_and_real = euclidean(result.betas, [20.0, 10.0])

        @test result isa EiveResult
        @test dist_clean_and_real < dist_dirty_and_real
        @test dist_eive_and_real < dist_dirty_and_real

    end


end
