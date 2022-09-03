using Test
using ErrorsInVariables
import Random

function euclidean(v1::Vector, v2::Vector)
    (v1 .- v2) .^ 2.0 |> sum |> sqrt
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

    result = eive(dirtyx = dirtyx, y = y, otherx = nothing, popsize = 100, numdummies = 10)

    dist_clean_and_real = euclidean(cleanbetas, [20.0, 10.0])
    dist_dirty_and_real = euclidean(dirtybetas, [20.0, 10.0])
    dist_eive_and_real = euclidean(result.betas, [20.0, 10.0])

    @test result isa EiveResult
    @test dist_clean_and_real < dist_dirty_and_real
    @test dist_eive_and_real < dist_dirty_and_real

end




@testset "Multiple regression (2 variables)" begin

    rng = Random.MersenneTwister(1234)

    n = 30

    deltax = randn(rng, n) * sqrt(3.0)

    cleanx = randn(rng, n) * sqrt(7.0)
    cleanx2 = randn(rng, n) * sqrt(7.0)

    e = randn(rng, n) * sqrt(5.0)

    y = 20.0 .+ 10.0 .* cleanx .+ 15.0 .* cleanx2 .+ e

    dirtyx = cleanx .+ deltax

    Xd = hcat(ones(n), dirtyx, cleanx2)
    Xc = hcat(ones(n), cleanx, cleanx2)
    dirtybetas = Xd \ y
    cleanbetas = Xc \ y

    result = eive(dirtyx = dirtyx, y = y, otherx = cleanx2, popsize = 100, numdummies = 10)

    dist_clean_and_real = euclidean(cleanbetas, [20.0, 10.0, 15.0])
    dist_dirty_and_real = euclidean(dirtybetas, [20.0, 10.0, 15.0])
    dist_eive_and_real = euclidean(result.betas, [20.0, 10.0, 15.0])

    @test result isa EiveResult
    @test dist_clean_and_real < dist_dirty_and_real
    @test dist_eive_and_real < dist_dirty_and_real

end




@testset "Multiple regression (3 variables)" begin

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

    result = eive(
        dirtyx = dirtyx,
        y = y,
        otherx = hcat(cleanx2, cleanx3),
        popsize = 100,
        numdummies = 10,
    )

    dist_clean_and_real = euclidean(cleanbetas, [20.0, 10.0, 15.0, 13.0])
    dist_dirty_and_real = euclidean(dirtybetas, [20.0, 10.0, 15.0, 13.0])
    dist_eive_and_real = euclidean(result.betas, [20.0, 10.0, 15.0, 13.0])

    @test result isa EiveResult
    @test dist_clean_and_real < dist_dirty_and_real
    @test dist_eive_and_real < dist_dirty_and_real

    #@info cleanbetas
    #@info dirtybetas
    #@info result
end
