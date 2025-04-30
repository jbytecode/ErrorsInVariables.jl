using Test 
using ErrorsInVariables
using Random


@testset "DD estimator" verbose=true begin 

    @testset "Simple Example" verbose=true begin 

        eps = 0.0001

        rng = Random.MersenneTwister(1234)

        n = 30

        deltax = randn(rng, n) * sqrt(3.0)

        cleanx = randn(rng, n) * sqrt(7.0)
        cleanx2 = randn(rng, n) * sqrt(7.0)
        cleanx3 = randn(rng, n) * sqrt(7.0)

        e = randn(rng, n) * sqrt(5.0)

        y = 20.0 .+ 10.0 .* cleanx .+ 15.0 .* cleanx2 .+ 13.0 .* cleanx3 .+ e

        dirtyx = cleanx .+ deltax

        Xd = hcat(ones(n), dirtyx, cleanx2, cleanx3)
        Xc = hcat(ones(n), cleanx, cleanx2, cleanx3)

        dirtybetas = Xd \ y

        cleanbetas = Xc \ y

        Z = dd(hcat(dirtyx, cleanx2, cleanx3), y)

        result = iv(Xd, y, Z)

        expected_dd_betas = [18.488555621677015, 6.283197918953448, 15.765140901269259, 12.061208813549317]

        for i in eachindex(expected_dd_betas)
            @test isapprox(result[i], expected_dd_betas[i], atol=eps)
        end

    end 
end 