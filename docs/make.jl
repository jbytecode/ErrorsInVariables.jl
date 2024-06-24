using Documenter, ErrorsInVariables

makedocs(
         format = Documenter.HTML(
                                  prettyurls = get(ENV, "CI", nothing) == "true",
                                  collapselevel = 2,
                                  # assets = ["assets/favicon.ico", "assets/extra_styles.css"],
                                 ),
         sitename="Errors-in-variables estimation using Compact Genetic Algorithms in Julia",
         authors = "Mehmet Hakan Satman",
         pages = [
                  "The Estimator" => "estimator.md",
                  "Multivariate Case" => "multivariate.md",
                  "CGA" => "cga.md", 
                  "Orthogonal Regression" => "orthogonalregression.md",
                 ]
        )


deploydocs(
    repo = "github.com/jbytecode/ErrorsInVariables.jl",
)