module Eive

include("cga.jl")
include("estimator.jl")



export CGA
export Estimator

import .Estimator: eive, EiveResult
import .CGA: cga, cgasample



export eive
export EiveResult
export cga
export cgasample


end # module
