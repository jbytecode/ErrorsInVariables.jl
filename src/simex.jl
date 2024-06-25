module SimulationExtrapolation

using Optim 


import ..SimpleEiveResult

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
    errvariance::Float64; numsims::Int = 1000)::SimpleEiveResult

    n, p = size(X)

    betamatrix = simex_multiple_iterations(X, y, λ, errorvarindex, 
                                                errvariance, numsims=numsims)

    # Extrapolate for λ = -1
    correctedbeta = extrapolate(λ, betamatrix, errorvarindex)
    
    # Remove the effect of the error variable from y vector 
    newy = y - X[:, errorvarindex] * correctedbeta
    
    # Regress newy to X matrix but without the error variable
    otherbetas = X[:, setdiff(1:p, errorvarindex)] \ newy

    # Combine correctedbeta and otherbetas
    result = zeros(p)
    result[errorvarindex] = correctedbeta
    result[setdiff(1:p, errorvarindex)] = otherbetas
     
    return SimpleEiveResult(result, true)
end 

function simex end 


end # end of module