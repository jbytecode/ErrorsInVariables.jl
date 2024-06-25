module SimulationExtrapolation

using Optim 

export simex 
export simex_single_iteration
export simex_multiple_iterations
export extrapolate


function simex_single_iteration(X::Matrix, 
                                y::Vector, 
                                λ::Float64, 
                                errorvarindex::Int, 
                                errvariance::Float64; numsims::Int = 1000)

    n, p = size(X)

    if errorvarindex > p
        throw(ArgumentError("errorvarindex must be less than or equal to the number of columns in X"))
    end

    if errvariance < 0
        throw(ArgumentError("errvariance must be non-negative"))
    end

    # Allocate memory for the results
    betasmatrix = zeros(numsims, p)

    XTemp = copy(X)

    part = sqrt(λ) * sqrt(errvariance)

    for current_iter in 1:numsims
        # Generate the error vector
        error =  part .* randn(n)

        # Add the error to the response variable
        XTemp[:, errorvarindex] .= X[:, errorvarindex] .+ error    

        # Fit the model
        β = XTemp \ y

        # Store the results
        betasmatrix[current_iter, :] = β
    end

    # Calculate column means of betasmatrix 
    betameans = sum(betasmatrix, dims=1) / numsims

    return betameans
end 



function simex_multiple_iterations(X::Matrix, 
    y::Vector, 
    λ::Vector, 
    errorvarindex::Int, 
    errvariance::Float64; numsims::Int = 1000)

    n, p = size(X)

    nlambdas = length(λ)

    lambdaresults = zeros(nlambdas, p)


    for i in 1:nlambdas
        lambdaresults[i, :] = simex_single_iteration(X, y, λ[i], errorvarindex, errvariance, numsims = numsims)
    end

    return lambdaresults
end 



function extrapolate(λ::Vector, betas::Matrix, errorvarindex::Int)::Float64

    x = λ
    
    y = betas[:, errorvarindex]

    function objective(par::Vector)::Float64
        yhat = par[1] .+ (par[2] ./ (par[3] .+ x))
        return sum((y .- yhat) .^ 2)
    end 

    initialpoint = rand(3)

    optresult = Optim.optimize(objective, initialpoint, Optim.Newton())

    par = Optim.minimizer(optresult)

    extrapol = par[1] + (par[2] / (par[3] + (-1)))

    return extrapol
end 


function simex(X::Matrix, 
    y::Vector, 
    λ::Vector, 
    errorvarindex::Int, 
    errvariance::Float64; numsims::Int = 1000)

    betamatrix = simex_multiple_iterations(x, y, lambdas, errorvarindex, 
                                                errvariance, numsims=numsims)

    # Extrapolate for λ = -1
    result = extrapolate(lambdas, betamatrix, errorvarindex) 

    return result
end 

function simex end 


end # end of module