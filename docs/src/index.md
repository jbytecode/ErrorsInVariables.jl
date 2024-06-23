[![Doc](https://img.shields.io/badge/docs-dev-blue.svg)](https://jbytecode.github.io/ErrorsInVariables.jl/)
[![codecov](https://codecov.io/gh/jbytecode/eive.jl/branch/main/graph/badge.svg?token=KMF7H1DS01)](https://codecov.io/gh/jbytecode/eive.jl)

# ErrorsInVariables.jl
Error-in-variables estimation using Compact Genetic Algorithms in Julia


# Usage 

For the single variable case 

```Julia 
julia> eive(dirtyx = dirtyx, y = y, otherx = nothing) 
```

and for the multiple regression 

```Julia 
julia> eive(dirtyx = dirtyx, y = y, otherx = matrixofotherx) 
```

Note that the method assumes there is only one erroneous variable in the set of independent variables.

# Example 

Lets generate data from the model $y = 20 + 10x^* + \varepsilon$

```julia
import Random
using ErrorsInVariables

rng = Random.MersenneTwister(1234)
n = 30
deltax = randn(rng, n) * sqrt(3.0)
cleanx = randn(rng, n) * sqrt(7.0)
e = randn(rng, n) * sqrt(5.0)
y = 20.0 .+ 10.0 .* cleanx .+ e
dirtyx = cleanx
eive(dirtyx = dirtyx, y = y, otherx = nothing) 
```

The result is 

```Julia
EiveResult([20.28458307772922, 9.456757289676714])
```

whereas OLS estimates are

```Julia
julia> X = hcat(ones(n), dirtyx);

julia> X \ y
2-element Vector{Float64}:
 17.94867860059858
  5.8099584879737876
```

and clearly biased towards to zero.

# Citation

```latex
@article{satman2015reducing,
  title={Reducing errors-in-variables bias in linear regression using compact genetic algorithms},
  author={Satman, M Hakan and Diyarbakirlioglu, Erkin},
  journal={Journal of Statistical Computation and Simulation},
  volume={85},
  number={16},
  pages={3216--3235},
  year={2015},
  publisher={Taylor \& Francis}
}

@article{Satman_Diyarbakirlioglu_2024, 
title={A Solution to Errors-in-variables Bias in Multivariate Linear Regression using Compact Genetic Algorithms},
volume={4}, 
url={https://journals.gen.tr/index.php/jame/article/view/2293}, 
DOI={10.53753/jame.2293}, 
number={1}, 
journal={JOURNAL OF APPLIED MICROECONOMETRICS}, 
author={Satman, Mehmet Hakan and Diyarbakırlıoğlu, Erkin}, 
year={2024}, month={Jun.}, 
pages={31–64} 
}
```