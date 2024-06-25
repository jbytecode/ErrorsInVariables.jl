module ErrorsInVariables

include("cga.jl")
include("estimator.jl")
include("meive.jl")
include("orthogonalregression.jl")
include("deming.jl")
include("simex.jl")


export CGA
export Estimator
export SimulationExtrapolation

import .Estimator: eive, EiveResult
import .Eivem: meive
import .CGA: cga, cgasample, converged 
import .OrthogonalRegression: orthogonal_regression
import .DemingRegression: deming
import .SimulationExtrapolation: simex


export eive
export meive
export EiveResult
export cga
export cgasample
export converged 
export orthogonal_regression
export deming 
export simex

end # module
