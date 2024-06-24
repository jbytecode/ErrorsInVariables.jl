using Test
using ErrorsInVariables
import Random


@testset "Multivariate CGA estimator" verbose = true begin 

@testset "Multivariate case - 2 response" begin
    rng = Random.MersenneTwister(1234)

    n = 30

    deltax = randn(rng, n) * sqrt(3.0)

    cleanx = randn(rng, n) * sqrt(7.0)

    e1 = randn(rng, n) * sqrt(5.0)
    e2 = randn(rng, n) * sqrt(5.0)

    y1 = 20.0 .+ 10.0 .* cleanx .+ e1
    y2 = 10.0 .+ 20.0 .* cleanx .+ e2

    dirtyx = cleanx .+ deltax

    Xd = hcat(ones(n), dirtyx)
    Xc = hcat(ones(n), cleanx)
    dirtybetas = Xd \ y1
    cleanbetas = Xc \ y1

    result = meive(
        dirtyx = dirtyx,
        y = hcat(y1, y2),
        otherx = nothing,
        popsize = 100,
        numdummies = 10,
    )

    dist_clean_and_real = euclidean(cleanbetas, [20.0, 10.0])
    dist_dirty_and_real = euclidean(dirtybetas, [20.0, 10.0])
    dist_eive_and_real = euclidean(result.betas, [20.0, 10.0])

    @test result isa EiveResult
    @test dist_clean_and_real < dist_dirty_and_real
    @test dist_eive_and_real < dist_dirty_and_real
end


end 