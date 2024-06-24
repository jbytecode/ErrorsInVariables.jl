using Test
using ErrorsInVariables
using Random 

@testset "Deming Regression" verbose=true begin 

@testset "R data example" verbose = true begin

    eps = 0.0001

    x = Float64[7, 8.3, 10.5, 9, 5.1, 8.2, 10.2, 10.3, 7.1, 5.9]
    y = Float64[7.9, 8.2, 9.6, 9, 6.5, 7.3, 10.2, 10.6, 6.3, 5.2]


    # The lambda parameter is set to 4.0
    # Lambda is the ratio of the variance of the errors in the y variable 
    # to the variance of the errors in the x variable.
    λ = 1.0 / 4.0

    result = deming(x, y, λ)

    @test result.converged

    @test isapprox(result.betas[1], -0.08974, atol = eps)
    @test isapprox(result.betas[2], 1.00119, atol = eps)
end



@testset "Simple Model" verbose = true begin

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

    λ = 1.0
    result = deming(dirtyx, y, λ)

    dist_clean_and_real = euclidean(cleanbetas, [20.0, 10.0])
    dist_dirty_and_real = euclidean(dirtybetas, [20.0, 10.0])
    dist_eive_and_real = euclidean(result.betas, [20.0, 10.0])


    @test result isa EiveResult
    @test dist_clean_and_real < dist_dirty_and_real
    @test dist_eive_and_real < dist_dirty_and_real

end


end 