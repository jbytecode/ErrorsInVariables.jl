module ErrorsInVariables

include("cga.jl")
include("estimator.jl")
include("meive.jl")
include("orthogonalregression.jl")
include("deming.jl")


export CGA
export Estimator


import .Estimator: eive, EiveResult
import .Eivem: meive
import .CGA: cga, cgasample, converged 
import .OrthogonalRegression: orthogonal_regression
import .DemingRegression: deming


export eive
export meive
export EiveResult
export cga
export cgasample
export converged 
export orthogonal_regression
export deming 


end # module
