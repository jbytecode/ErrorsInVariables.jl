module DemingRegression 

import ..Estimator: EiveResult

export deming 


"""
    deming(x::Vector, y::Vector, λ::Float64 = 1.0)::EiveResult

# Description

    This function estimates the parameters of a linear model using the Deming regression method.

# Arguments

- `x::Vector`: A vector of the independent variable.
- `y::Vector`: A vector of the dependent variable.
- `λ::Float64 = 1.0`: The ratio of the variance of the errors in the y variable to the variance of the errors in the x variable.

# Returns

- `EiveResult`: An object that contains the estimated parameters of the model.

# Example
```julia
using ErrorsInVariables

x = Float64[7, 8.3, 10.5, 9, 5.1, 8.2, 10.2, 10.3, 7.1, 5.9]
y = Float64[7.9, 8.2, 9.6, 9, 6.5, 7.3, 10.2, 10.6, 6.3, 5.2]

# The lambda parameter is set to 4.0
# Lambda is the ratio of the variance of the errors in the y variable 
# to the variance of the errors in the x variable.
λ = 1.0 / 4.0

result = deming(x, y, λ)
``` 

# References
- https://en.wikipedia.org/wiki/Deming_regression
"""
function deming(x::Vector, y::Vector, λ::Float64 = 1.0)::EiveResult
   
    n = length(x)

    if length(y) != n
        throw(ArgumentError("X and y must have the same length"))
    end

    xbar = sum(x) / n 
    ybar = sum(y) / n


    xmxbar = (x .- xbar)
    
    ymybar = (y .- ybar)

    sxx = sum(xmxbar .^ 2) / n
    
    syy = sum(ymybar .^ 2) / n 

    sxy = sum(xmxbar .* ymybar) / n

    b1 = (syy - λ * sxx + sqrt((syy - λ * sxx) ^ 2 + 4 * λ * sxy ^ 2)) / (2 * sxy)

    b0 = ybar - b1 * xbar

    return EiveResult([b0, b1], true)
end # end of function


end # end of module