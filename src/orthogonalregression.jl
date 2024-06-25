module OrthogonalRegression

using Optim

export orthogonal_regression

import ..SimpleEiveResult

"""
    orthogonal_regression(
        X::Matrix,
        y::Vector;
        xhasintercept::Bool=true,
        maxiterations::Int=10000,
        initialbetas::Union{Nothing, Vector} = nothing)::SimpleEiveResult

# Description:

The function performs orthogonal regression. The method searches for a set of 
parameters that minimize the sum of squares of orthogonal residuals.

# Arguments:

`X::Matrix`: Independent variables. The first column should be 1 for the intercept.

`y::Vector`: Dependent variable.

`xhasintercept::Bool`: If true, the first column of X is considered as the intercept.

`maxiterations::Int`: Maximum number of iterations.

`initialbetas::Union{Nothing, Vector}`: Initial values for the parameters. 
If nothing, random values are used.


# Returns:

`::EiveResult`: A struct containing the estimated parameters and convergence status.

# Examples

```julia
X = Float64[1 1.0; 1 0.6; 1 1.2; 1 1.4; 1 0.2]
y = Float64[0.5, 0.3, 0.7, 1.0, 0.2]

result = orthogonal_regression(X, y)
```

# References
- https://davegiles.blogspot.com/2014/11/orthogonal-regression-first-steps.html
"""
function orthogonal_regression(X::Matrix,
    y::Vector;
    xhasintercept::Bool=true,
    maxiterations::Int=10000,
    initialbetas::Union{Nothing, Vector} = nothing)::SimpleEiveResult

    _, p = size(X)

    if isnothing(initialbetas)
        initial = rand(p)
    else
        initial = initialbetas
    end

    function objective(params)
        res = y .- X * params
        exclude_intercept = params[2:end]
        if xhasintercept
            exclude_intercept = params[2:end]
            ortres = res ./ sqrt(1 .+ sum(exclude_intercept .^ 2))
        else
            ortres = res ./ sqrt(1 .+ sum(params .^ 2))
        end
        retval = sum(ortres .^ 2)
        return retval
    end

    result = Optim.optimize(objective,
        initial,
        Newton(),
        Optim.Options(g_tol=1e-6, iterations=maxiterations))

    return SimpleEiveResult(
        Optim.minimizer(result), 
        Optim.converged(result))
end

end # end of module