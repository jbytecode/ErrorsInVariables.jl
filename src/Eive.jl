module Eive

include("cga.jl")
include("estimator.jl")



export CGA 
export Estimator 

import .Estimator: eive 
import .CGA: cga, cgasample 



export eive
export cga
export cgasample


end # module
