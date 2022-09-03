module Eives

include("cga.jl")
include("estimator.jl")



export CGA
export Estimator

import .Estimator: eives, EivesResult
import .CGA: cga, cgasample



export eives
export EivesResult
export cga
export cgasample


end # module
