module ErrorsInVariables

include("cga.jl")
include("estimator.jl")
include("meive.jl")


export CGA
export Estimator


import .Estimator: eive, EiveResult
import .Eivem: meive
import .CGA: cga, cgasample, converged 



export eive
export meive
export EiveResult
export cga
export cgasample
export converged 


end # module
