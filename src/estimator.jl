module Estimator

import Eive.CGA: cga

function eive(
    dirtyx::Array{T,1},
    y::Array{T,1},
    otherx::Union{Nothing, Array{T,2}, Array{T, 1}},
    popsize::Int = 50,
    numdummies::Int = 10,
) where {T<:Real}

    if isnothing(otherx)
        eivewithoutotherx(dirtyx, y, popsize, numdummies)
    else
        return eivewithotherx(dirtyx, y, otherx, popsize, numdummies)
    end
end

function eivewithotherx(
    dirtyx::Array{T,1},
    y::Array{T,1},
    otherx::Union{Array{T,2}, Array{T, 1}},
    popsize::Int = 50,
    numdummies::Int = 10,
) where {T<:Real}


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

    finalbits = cga(chsize = chsize, costfunction = costfn, popsize = popsize)

    auxX = reshape(finalbits, n, numdummies)
    betas = auxX \ dirtyx
    cleanX = auxX * betas

    X = hcat(myones, cleanX, otherx)

    outerbetas = X \ y

    return outerbetas
end


function eivewithoutotherx(
    dirtyx::Array{T,1},
    y::Array{T,1},
    popsize::Int = 50,
    numdummies::Int = 10,
) where {T<:Real}


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

    finalbits = cga(chsize = chsize, costfunction = costfn, popsize = popsize)

    auxX = reshape(finalbits, n, numdummies)
    betas = auxX \ dirtyx
    cleanX = auxX * betas

    X = hcat(myones, cleanX)

    outerbetas = X \ y

    return outerbetas
end


export eive

end # end of module
