module ErrorsInVariables

using SumTypes

@sum_type Either{A, B} begin
    Left{A}(::A)
    Right{B}(::B)
end

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

export Either, Left, Right

export CGA
export Estimator
export SimulationExtrapolation
export IV
export DD

import .SumTypes: @cases, @sum_type
import .Estimator: eive
import .Eivem: meive
import .CGA: cga, cgasample, converged 
import .OrthogonalRegression: orthogonal_regression
import .DemingRegression: deming
import .SimulationExtrapolation: simex, simex_single_iteration, simex_multiple_iterations, extrapolate
import .IV: iv
import .DD: dd 

export @cases, @sum_type

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
