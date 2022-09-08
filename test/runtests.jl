using Test
using ErrorsInVariables


function euclidean(v1::Vector, v2::Vector)
    (v1 .- v2) .^ 2.0 |> sum |> sqrt
end

include("testcga.jl")
include("testestimator.jl")
include("testmeive.jl")
