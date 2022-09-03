using Test
using Eive


@testset "sample" begin
    @testset "sample full of ones" begin
        @test cgasample([1, 1, 1, 1, 1]) == [1, 1, 1, 1, 1]
    end

    @testset "sample full of zeros" begin
        @test cgasample([0, 0, 0, 0, 0]) == [0, 0, 0, 0, 0]
    end

    @testset "sample with half ones" begin
        smp = cgasample([0.5, 0.5, 0.5, 0.5, 0.5])
        for element in smp
            @test (element >= 0.0) || (element <= 1.0)
        end
    end
end

@testset "cga" begin
    function f(bits)
        return sum(bits)
    end

    result = cga(chsize = 10, costfunction = f, popsize = 100)

    for element in result
        @test element == 0
    end
end
