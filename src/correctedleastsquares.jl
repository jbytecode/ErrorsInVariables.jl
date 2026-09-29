module CorrectedLeastSquares

export corrected_least_squares

import LinearAlgebra: Symmetric, eigvals, opnorm, rank
import ..SimpleEiveResult

"""
    corrected_least_squares(
        X::Matrix,
        y::Vector,
        errorcovariance::Matrix,
    )::SimpleEiveResult

Estimate linear-regression coefficients with the corrected least-squares (CLS)
method when one or more columns of `X` contain additive measurement error.

## Definition

Let the observed design matrix satisfy ``X = X^* + \\Delta``, where ``X^*`` is
the unobserved error-free design matrix and each row of ``\\Delta`` has known
covariance matrix ``\\Sigma_\\Delta``. Under independent, mean-zero measurement
errors, CLS corrects the biased normal equations as

```math
\\hat{\\beta}_{CLS} =
    \\left(X^\\top X - n\\Sigma_\\Delta\\right)^{-1}X^\\top y.
```

`errorcovariance` is ``\\Sigma_\\Delta`` and must have one row and column per
column of `X`. If `X` includes an intercept column, its corresponding row and
column in `errorcovariance` must be zero.

## Arguments

- `X`: Observed `n × p` design matrix.
- `y`: Response vector with `n` observations.
- `errorcovariance`: Known or externally estimated `p × p` covariance matrix
  of the row-wise measurement errors in `X`.

## Returns

A `SimpleEiveResult` containing the corrected coefficient vector. CLS is a
closed-form estimator, so `converged` is always `true`.

## References

- Fuller, W. A. (1987). *Measurement Error Models*. Wiley, Sections 1.3--1.4.
- Carroll, R. J., Ruppert, D., Stefanski, L. A., & Crainiceanu, C. M. (2006).
  *Measurement Error in Nonlinear Models: A Modern Perspective* (2nd ed.).
  Chapman & Hall/CRC, Chapter 3.
"""
function corrected_least_squares(
    X::Matrix,
    y::Vector,
    errorcovariance::Matrix,
)::SimpleEiveResult
    n, p = size(X)

    n > 0 || throw(ArgumentError("X must contain at least one observation"))
    length(y) == n || throw(ArgumentError("X and y must have the same number of observations"))
    size(errorcovariance) == (p, p) ||
        throw(ArgumentError("errorcovariance must have size ($p, $p)"))
    all(isfinite, X) || throw(ArgumentError("X must contain only finite values"))
    all(isfinite, y) || throw(ArgumentError("y must contain only finite values"))
    all(isfinite, errorcovariance) ||
        throw(ArgumentError("errorcovariance must contain only finite values"))

    tolerance = sqrt(eps(Float64)) * max(1.0, opnorm(errorcovariance, Inf))
    isapprox(errorcovariance, transpose(errorcovariance); atol=tolerance, rtol=0.0) ||
        throw(ArgumentError("errorcovariance must be symmetric"))
    minimum(eigvals(Symmetric(errorcovariance))) >= -tolerance ||
        throw(ArgumentError("errorcovariance must be positive semidefinite"))

    corrected_crossproduct = X' * X - n * errorcovariance
    rank(corrected_crossproduct) == p ||
        throw(ArgumentError("corrected normal-equation matrix is rank deficient"))

    return SimpleEiveResult(
        Vector{Float64}(corrected_crossproduct \ (X' * y)),
        true,
    )
end

end # module
