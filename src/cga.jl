module CGA

import Random: AbstractRNG, MersenneTwister

"""
Generates a binary vector of values using a probability vector.
Each single element of the probability vector is the probability of bit having 
the value of 1. When the probability vector is [1, 1, 1, ..., 1] then the sampled
vector is [1.0, 1.0, 1.0, ..., 1.0] whereas it is [0.0, 0.0, 0.0, ..., 0.0] when the probability vector
is a vector of zeros. The CGA (compact genetic algorithms) search is started using the 
probability vector of [0.5, 0.5, 0.5, ..., 0.5] which produces random vectors of either
zeros or ones.
 
# Examples
```jldoctest
julia> sample([1, 1, 1, 1, 1])
5-element Vector{Int}:
 1
 1
 1
 1
 1
julia> cgasample(ones(10) * 0.5)
10-element Vector{Int}:
 1
 1
 1
 1
 0
 0
 0
 1
 1
 0
```
"""
function cgasample(
	probvector::Vector{T}, 
	rng::RNGType)::Vector{Int} where {T <: Real, RNGType <: AbstractRNG}

	function sampler(x)
		if rand(rng) < x
			1
		else
			0
		end
	end

	return map(sampler, probvector)
end

function cgasampleinto!(
	probvector::Vector{T}, 
	rng::RNGType, 
	samplevector::Vector{Int})::Nothing where {T <: Real, RNGType <: AbstractRNG}

	for i ∈ 1:length(probvector)
		if rand(rng) < probvector[i]
			samplevector[i] = 1
		else
			samplevector[i] = 0
		end
	end

	return nothing
end



function converged(probvector::Vector{T}) where {T <: Real}
	return all(x -> iszero(x) || isone(x), probvector)
end


function cgaupdate!(
	probvector::Vector{T},
	winner::Vector{Int},
	loser::Vector{Int},
	mutation::Float64)::Nothing where {T <: Real}

	chsize = length(probvector)

	for i ∈ 1:chsize
		if winner[i] != loser[i]
			if isone(winner[i])
				probvector[i] = min(1.0, probvector[i] + mutation)
			else
				probvector[i] = max(0.0, probvector[i] - mutation)
			end
		end
	end
end

"""
Performs a CGA (Compact Genetic Algorithm) search for minimization of an objective function.
In the example below, the objective function is to minimize sum of bits of a binary vector.
The search method results the optimum vector of [0, 0, ..., 0] where the objective function is zero.


# Examples
```jldoctest
julia> function f(x)
		   return sum(x)
	   end
f (generic function with 1 method)
julia> cga(chsize = 10, costfunction = f, popsize = 100)
10-element Vector{Int}:
 0
 0
 0
 0
 0
 0
 0
 0
 0
 0
```
"""
function cga(;
	chsize::Int,
	costfunction::Function,
	popsize::Int,
	rng::RNGType = MersenneTwister(0),
)::Vector{Int} where {RNGType<:AbstractRNG}

	probvector = ones(Float64, chsize) * 0.5
	mutation = 1.0 / convert(Float64, popsize)
	candidate1 = Array{Int, 1}(undef, chsize)
	candidate2 = Array{Int, 1}(undef, chsize)

	while !converged(probvector)

		cgasampleinto!(probvector, rng, candidate1)
		cgasampleinto!(probvector, rng, candidate2)

		cost1 = costfunction(candidate1)
		cost2 = costfunction(candidate2)

		if cost1 < cost2 
			cgaupdate!(probvector, candidate1, candidate2, mutation)
		else
			cgaupdate!(probvector, candidate2, candidate1, mutation)
		end

	end
	
	cgasampleinto!(probvector, rng, candidate1)
	return candidate1
end

end # End of module
