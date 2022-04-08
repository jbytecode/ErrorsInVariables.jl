module Eive

include("cga.jl")
include("estimator.jl")



export CGA 
export Estimator 

import .Estimator: eive 

export eive



end # module
