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
function dd(xvec::Vector, yvec::Vector)::Matrix 

    n = length(xvec)

    IK = ones(n)

    x = xvec .- mean(xvec)

    y = yvec .- mean(yvec)

    z1 = x .* x 

    z2 = x .* y 

    z3 = y .* y

    z4 = x .* x .* x .- 3 .* x .* (mean(x'x / n) .* IK)

    # z5 = x .* x .* y .- 2 .* x .* (mean(y'x / n) .* IK) .- y * (IK' * mean(x'x / n) .* IK)

    return hcat(z1, z2, z3, z4)

end 

end # End of module DD