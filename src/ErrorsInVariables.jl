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


export CGA
export Estimator
export SimulationExtrapolation

import .Estimator: eive
import .Eivem: meive
import .CGA: cga, cgasample, converged 
import .OrthogonalRegression: orthogonal_regression
import .DemingRegression: deming
import .SimulationExtrapolation: simex


export eive
export meive
export EiveResult, SimpleEiveResult
export cga
export cgasample
export converged 
export orthogonal_regression
export deming 
export simex

end # module
