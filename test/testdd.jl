using Test 
using ErrorsInVariables


@testset "DD estimator" verbose=true begin 

    @testset "Simple Example" verbose=true begin 

        X = [1.0, 2.0, 3.0, 4.0, 5.0, 6, 7, 8]

        y = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0]

        result = dd(X, y)

        println("-----------------------")
        display(result)
    end 
end 