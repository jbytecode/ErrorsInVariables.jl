using Test, Random 
using ErrorsInVariables
using ErrorsInVariables.SimulationExtrapolation

@testset "Simex" verbose = true begin

    @testset "Simex Single Iteration for λ = 0" begin

        eps = 0.0001

        x = Float64[1 1; 1 2; 1 3; 1 4; 1 5]
        y = Float64[1, 2, 3, 4, 5]

        λ = 0.0

        errorvarindex = 2

        errvariance = 1.0

        betameans = simex_single_iteration(x, y, λ, errorvarindex, errvariance, numsims=1000)

        olsresult = x \ y

        @test isapprox(betameans[1], olsresult[1], atol=eps)
        @test isapprox(betameans[2], olsresult[2], atol=eps)

    end



    @testset "Simex Single Iteration for λ = 0.5" begin

        x = Float64[1 1; 1 2; 1 3; 1 4; 1 5]
        y = Float64[1, 2, 3, 4, 5]

        λ = 0.1

        errorvarindex = 2

        errvariance = 1.0

        betameans = simex_single_iteration(x, y, λ, errorvarindex, errvariance, numsims=1000)

        olsresult = x \ y

        @test betameans[1] != olsresult[1]
        @test betameans[2] != olsresult[2]

    end


    @testset "Simex Multiple iterations for several λ values" begin

        eps = 0.0001


        x = Float64[1 1; 1 2; 1 3; 1 4; 1 5]
        y = Float64[2, 4, 6, 8, 10]


        errorvarindex = 2

        errvariance = 1.0

        lambdas = [i / 8 for i in 0:16]

        result = simex_multiple_iterations(x, y, lambdas, errorvarindex, errvariance, numsims=1000)

        olsresult = x \ y

        @test isapprox(result[1, 1], olsresult[1], atol=eps)
        @test isapprox(result[1, 2], olsresult[2], atol=eps)

        @test result isa Matrix

        @test size(result) == (length(lambdas), 2)
    end


    @testset "Simex Extrapolate for λ = -1" begin

        eps = 0.1

        numsims = 100000

        x = Float64[1 1; 1 2; 1 3; 1 4; 1 5]
        y = Float64[2, 4, 6, 8, 10]

        errorvarindex = 2

        errvariance = 1.0

        lambdas = [i / 8 for i in 0:16]

        betamatrix = simex_multiple_iterations(x, y, lambdas, errorvarindex, 
                                                errvariance, numsims=numsims)

        # Extrapolate for λ = -1
        result = extrapolate(lambdas, betamatrix, errorvarindex)

        @test isapprox(result, 2.5, atol=eps)

    end


    @testset "Simple Model Extrapolation" verbose = true begin

        eps = 0.0001

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
    
        errorvarindex = 2

        # We assume that the error in the x variable is 3.0
        # Because of the data generating process is known here.
        errvariance = 3.0

        lambdas = [i / 8 for i in 0:16]

        betamatrix = simex_multiple_iterations(Xd, y, lambdas, errorvarindex, 
                                                errvariance, numsims=10000)

        # Extrapolate for λ = -1
        result = extrapolate(lambdas, betamatrix, errorvarindex)

        @test isapprox(result, 8, atol=1)
     
    end
end