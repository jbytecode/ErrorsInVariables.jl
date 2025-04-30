module BitOptimizer


using Evolutionary
import ..CGA: cga

export cgaoptimizer
export gaoptimizer


function cgaoptimizer(; chsize::Int, costfunction::F, popsize::Int, rng) where {F<:Function}
    return cga(chsize=chsize, costfunction=costfunction, popsize=popsize, rng=rng)
end

function gaoptimizer(; chsize::Int, costfunction::F, popsize::Int, rng) where {F<:Function}

    optimizer = GA(populationSize = popsize,
                   crossoverRate = 0.9,
                   mutationRate = 0.5,
                   selection=tournament(2),
                   mutation=flip, 
                   crossover=UX,
                   epsilon = 2)
                   

    optimresult = Evolutionary.optimize(costfunction,
                            BitVector(zeros(chsize)),
                            optimizer,
                            Evolutionary.Options(iterations=10000))

    return (optimresult.minimizer)
end

end