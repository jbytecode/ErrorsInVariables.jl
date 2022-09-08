module Eivem


export meive

import ..Estimator: EiveResult
import ..CGA: cga
import Random: AbstractRNG, MersenneTwister


"""
    eive(;
    dirtyx::Array{T, 1},
    y::Array{T, 2},
    otherx::Union{Nothing,Array{T,2},Array{T,1}},
    popsize::Int = 50,
    numdummies::Int = 10,
    rng::AbstractRNG = MersenneTwister(1234)
)::EiveResult where {T<:Real}

# Description:
This is the multivariate case of eive(). Please see eive() function. 
In the multivariate case, the y is not vector, but a matrix of multiple 
or repeated measurements of the response variable. This can be considered 
as multivariate regression as well as regressions with repeated measurements.  

# Arguments:

- dirtyx: Independent variable measured with some error
- y: nxp matrix of dependent variables where n is the number of observations and p is the number of dependent variables
- otherx: Matrix of other independent variables
- popsize: Number of individuals in the population (optional)
- numdummies: Number of dummy variables to use (optional)
- rng: Random number generator (optional)

# Examples

```julia-repl
julia> import Random
julia> using ErrorsInVariables
julia> rng = Random.MersenneTwister(1234)
julia> n = 30
julia> deltax = randn(rng, n) * sqrt(3.0)
julia> cleanx = randn(rng, n) * sqrt(7.0)
julia> e1 = randn(rng, n) * sqrt(5.0)
julia> e2 = randn(rng, n) * sqrt(5.0)
julia> y1 = 20.0 .+ 10.0 .* cleanx .+ e1
julia> y2 = 10.0 .+ 15.0 .* cleanx .+ e2
julia> dirtyx = cleanx + deltax

julia> # Getting bias-reduced estimates
julia> meive(dirtyx = dirtyx, y = hcat(y1, y2), otherx = nothing) 

EiveResult([19.65449584842238, 9.21108792897651])

julia> X = hcat(ones(n), dirtyx);

julia> # Biased OLS estimates:
julia> X \\ y1
2-element Vector{Float64}:
17.94867860059858
  5.8099584879737876
```

# References 

Satman, M. Hakan, and Erkin Diyarbakirlioglu. "Reducing errors-in-variables bias in linear 
regression using compact genetic algorithms." Journal of Statistical Computation and Simulation
 85.16 (2015): 3216-3235.

"""
function meive(;
    dirtyx::Array{T, 1},
    y::Array{T, 2},
    otherx::Union{Nothing,Array{T,2},Array{T,1}},
    popsize::Int = 50,
    numdummies::Int = 10,
    rng::AbstractRNG = MersenneTwister(1234))::EiveResult where {T<:Real}

    if isnothing(otherx)
        return meivewithoutotherx(dirtyx, y, popsize, numdummies, rng)
    else
        return meivewithotherx(dirtyx, y, otherx, popsize, numdummies, rng)
    end
end


function meivewithotherx(dirtyx::Array{T,1},y::Array{T,2}, otherx::Union{Array{T,2},Array{T,1}}, popsize::Int = 50, numdummies::Int = 10, rng::AbstractRNG = MersenneTwister(1234))::EiveResult where {T<:Real}


    n = length(dirtyx)
    myones = ones(Float64, n)
    chsize = n * numdummies


    function costfn(bits::Array{Int,1})
        auxX = reshape(bits, n, numdummies)
        betas = auxX \ dirtyx
        cleanX = auxX * betas

        X = hcat(myones, cleanX, otherx)

        _ , np = size(y)
        totalres = 0.0
        for i in 1:np
            currenty = y[:, i]
            outerbetas = X \ currenty
            res = currenty .- X * outerbetas
            totalres += sum(res .^ 2.0) 
        end 
        return totalres
    end

    finalbits = cga(chsize = chsize, costfunction = costfn, popsize = popsize, rng = rng)

    auxX = reshape(finalbits, n, numdummies)
    betas = auxX \ dirtyx
    cleanX = auxX * betas

    X = hcat(myones, cleanX, otherx)

    outerbetas = X \ y[:,1]

    return EiveResult(outerbetas)
end


function meivewithoutotherx(dirtyx::Array{T,1},y::Array{T,2},popsize::Int = 50,numdummies::Int = 10,rng::AbstractRNG = MersenneTwister(1234))::EiveResult where {T<:Real}

    n = length(dirtyx)
    myones = ones(Float64, n)
    chsize = n * numdummies


    function costfn(bits::Array{Int,1})
        auxX = reshape(bits, n, numdummies)
        betas = auxX \ dirtyx
        cleanX = auxX * betas

        X = hcat(myones, cleanX)

        _ , np = size(y)
        totalres = 0.0
        for i in 1:np
            currenty = y[:, i]
            outerbetas = X \ currenty
            res = currenty .- X * outerbetas
            totalres += sum(res .^ 2.0)
        end 
        return totalres
    end

    finalbits = cga(chsize = chsize, costfunction = costfn, popsize = popsize, rng = rng)

    auxX = reshape(finalbits, n, numdummies)
    betas = auxX \ dirtyx
    cleanX = auxX * betas

    X = hcat(myones, cleanX)

    outerbetas = X \ y[:, 1]

    return EiveResult(outerbetas)
end


end # end of module eivem 