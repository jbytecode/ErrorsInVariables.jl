module ErrorsInVariables

include("cga.jl")
include("estimator.jl")
include("meive.jl")
include("orthogonalregression.jl")


export CGA
export Estimator


import .Estimator: eive, EiveResult
import .Eivem: meive
import .CGA: cga, cgasample, converged 
import .OrthogonalRegression: orthogonal_regression



export eive
export meive
export EiveResult
export cga
export cgasample
export converged 
export orthogonal_regression


end # module
