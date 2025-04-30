module IV

export iv

import ..ErrorsInVariables: Either, Left, Right


"""
    iv(X::Matrix, y::Vector, Z::Matrix)::Vector

Computes the Instrumental Variable (IV) estimator for a linear regression model.

# Description

If the number of instruments is equal to the number of predictors, the IV estimator is 
computed using the formula:
```julia
    betas = inv(Z'X)Z'y
```
If the number of instruments is greater than the number of predictors, the IV estimator is
computed using the formula:
```julia
    part = X' * Z * inv(Z' * Z) * Z' 
    betas = inv(part * X) * part * y
```
"""
function iv(X::Matrix, y::Vector, Z::Matrix)::Either

    numberofinstruments = size(Z, 2)
    
    numberofpredictors = size(X, 2)

    if numberofinstruments < numberofpredictors 
        return Left("Error: Number of instruments must be greater than or equal to the number of predictors.")
    end

    if numberofinstruments == numberofpredictors
        return Right(inv(Z'X)Z'y)
    else
        part = X' * Z * inv(Z' * Z) * Z' 
        betas = inv(part * X) * part * y
        return Right{Vector}(betas)
    end

end

end # end of module IV