using Test
using ErrorsInVariables

@testset "IV Estimator" verbose=true begin 

    @testset "OLS Test" verbose = true begin

        # Test case 1: Simple OLS regression
        # When the instrument matrix is the same as the predictor matrix
        # and the response variable is a linear function of the predictors
        # then the IV estimator should return the same coefficients as OLS.

        X = [1.0 1.0; 1.0 2.0; 1.0 3.0; 1.0 4.0; 1.0 5.0]
        Z = X 
        y = [2.0, 4.0, 6.0, 8.0, 10.0]

        expected_betas = [0.0, 2.0]

        ivbetas = iv(X, y, Z)

        @test ivbetas isa Vector

        @test length(ivbetas) == 2

        @test isapprox(ivbetas[1], expected_betas[1], atol=1e-5)

        @test isapprox(ivbetas[2], expected_betas[2], atol=1e-5)
    end

end

    