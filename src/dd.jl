module DD


export dd 

import LinearAlgebra: I
import Statistics: mean

"""

    dd(xvec::Vector, yvec::Vector)::Matrix

Computes the Higher Moment Estimator for a linear regression model with errors in the variables.

# Reference

- Dagenais, Marcel G., and Denyse L. Dagenais. "Higher moment estimators for linear regression models with errors in the variables." 
Journal of Econometrics 76.1-2 (1997): 193-221.
"""
function dd(xmat::Matrix, yvec::Vector)::Matrix 

    n, p = size(xmat)

    o = ones(p)
    IK = I(p)

    x = similar(xmat)
    for i in 1:p
        x[:, i] = xmat[:, i] .- mean(xmat[:, i])
    end

    y = yvec .- mean(yvec)

    z1 = x .* x 

    z2 = x .* y 

    z3 = y .* y

    z4 = x .* x .* x .- 3 .* x * (mean(x'x ./ n) * IK)

    z5 = x .* x .* y .- 2 .* x * (mean(y'x ./ n) .* IK) .- y .* (mean(x'x / n) .* ones(n))

    z6 = x .* y .* y .- x * (mean(y'y ./ n) .* IK) .- 2 * y .* (mean(y'x / n) .* ones(n))

    z7 = y .* y .* y .- 3 .* y .* (mean(y'y ./ n))

    return hcat(z1, z2, z3, z4, z5, z6, z7)

end 

end # End of module DD