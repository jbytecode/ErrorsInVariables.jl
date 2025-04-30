module ErrorsInVariables

abstract type EiveResult end 

struct SimpleEiveResult <: EiveResult
    betas::Vector{Float64}
    converged::Bool
end

include("cga.jl")
include("estimator.jl")
include("meive.jl")
include("orthogonalregression.jl")
include("deming.jl")
include("simex.jl")
include("iv.jl")
include("dd.jl")


export CGA
export Estimator
export SimulationExtrapolation
export IV
export DD


import .Estimator: eive
import .Eivem: meive
import .CGA: cga, cgasample, converged 
import .OrthogonalRegression: orthogonal_regression
import .DemingRegression: deming
import .SimulationExtrapolation: simex, simex_single_iteration, simex_multiple_iterations, extrapolate
import .IV: iv
import .DD: dd 


export eive
export meive

export EiveResult, SimpleEiveResult

export cga

export cgasample

export converged 

export orthogonal_regression

export deming 

export simex
export simex_single_iteration
export simex_multiple_iterations
export extrapolate

export iv 

export dd
end # module
