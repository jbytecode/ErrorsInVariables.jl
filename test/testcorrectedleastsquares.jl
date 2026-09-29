using Test
using ErrorsInVariables

@testset "Corrected least squares" verbose = true begin
    @testset "Recovers coefficients with known measurement error covariance" begin
        latent_x = Float64[-2, -1, 0, 1, 2]
        measurement_error = Float64[1, -2, 0, 2, -1]
        X = hcat(ones(5), latent_x .+ measurement_error)
        y = 3.0 .+ 4.0 .* latent_x
        errorcovariance = Float64[0 0; 0 2]

        result = corrected_least_squares(X, y, errorcovariance)

        @test result isa EiveResult
        @test result.converged
        @test result.betas ≈ [3.0, 4.0]
    end

    @testset "Matches OLS without measurement error" begin
        X = Float64[1 1; 1 2; 1 3; 1 4]
        y = Float64[3, 5, 7, 9]

        result = corrected_least_squares(X, y, zeros(2, 2))

        @test result.betas ≈ X \ y
    end

    @testset "Validates covariance matrix" begin
        X = Float64[1 1; 1 2; 1 3; 1 4]
        y = Float64[3, 5, 7, 9]

        @test_throws ArgumentError corrected_least_squares(X, y, zeros(3, 3))
        @test_throws ArgumentError corrected_least_squares(X, y, [0.0 1.0; 0.0 0.0])
        @test_throws ArgumentError corrected_least_squares(X, y, [0.0 0.0; 0.0 -1.0])
    end
end
