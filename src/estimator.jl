module Estimator

export eive
export EiveResult


import ..CGA: cga
import Random: AbstractRNG, MersenneTwister

struct EiveResult
    betas::Vector{Float64}
end

function eive(;
    dirtyx::Array{T,1},
    y::Array{T,1},
    otherx::Union{Nothing,Array{T,2},Array{T,1}},
    popsize::Int = 50,
    numdummies::Int = 10,
    rng::AbstractRNG = MersenneTwister(1234)
)::EiveResult where {T<:Real}

    if isnothing(otherx)
        eivewithoutotherx(dirtyx, y, popsize, numdummies, rng)
    else
        return eivewithotherx(dirtyx, y, otherx, popsize, numdummies, rng)
    end
end

function eivewithotherx(
    dirtyx::Array{T,1},
    y::Array{T,1},
    otherx::Union{Array{T,2},Array{T,1}},
    popsize::Int = 50,
    numdummies::Int = 10,
    rng::AbstractRNG = MersenneTwister(1234)
)::EiveResult where {T<:Real}


    n = length(dirtyx)
    myones = ones(Float64, n)
    chsize = n * numdummies


    function costfn(bits::Array{Int,1})
        auxX = reshape(bits, n, numdummies)
        betas = auxX \ dirtyx
        cleanX = auxX * betas

        X = hcat(myones, cleanX, otherx)

        outerbetas = X \ y
        res = y .- X * outerbetas
        return sum(res .^ 2.0)
    end

    finalbits = cga(chsize = chsize, costfunction = costfn, popsize = popsize, rng = rng)

    auxX = reshape(finalbits, n, numdummies)
    betas = auxX \ dirtyx
    cleanX = auxX * betas

    X = hcat(myones, cleanX, otherx)

    outerbetas = X \ y

    return EiveResult(outerbetas)
end


function eivewithoutotherx(
    dirtyx::Array{T,1},
    y::Array{T,1},
    popsize::Int = 50,
    numdummies::Int = 10,
    rng::AbstractRNG = MersenneTwister(1234)
)::EiveResult where {T<:Real}


    n = length(dirtyx)
    myones = ones(Float64, n)
    chsize = n * numdummies


    function costfn(bits::Array{Int,1})
        auxX = reshape(bits, n, numdummies)
        betas = auxX \ dirtyx
        cleanX = auxX * betas

        X = hcat(myones, cleanX)

        outerbetas = X \ y
        res = y .- X * outerbetas
        return sum(res .^ 2.0)
    end

    finalbits = cga(chsize = chsize, costfunction = costfn, popsize = popsize, rng = rng)

    auxX = reshape(finalbits, n, numdummies)
    betas = auxX \ dirtyx
    cleanX = auxX * betas

    X = hcat(myones, cleanX)

    outerbetas = X \ y

    return EiveResult(outerbetas)
end


end # end of module
