module Surveys
using Statistics, StatsBase, StatsAPI, StatsModels, DiffResults, ForwardDiff, PDMats
using LinearAlgebra, Distributions, GLM
include("docboilerplate.jl")

export SampleSum, SI, WithReplacement, WithoutReplacement, SurveyDesign

# Abstract Interfaces

"""
Abstract Survey Design. When passed as an argument to `sum` or `mean`, computes the design-based estimate and variance.
"""
abstract type SurveyDesign end

"""
With-replacement survey design. Used for probability-weighted ratio estimators (Hansen-Hurwitz) for population totals and means.
"""
struct WithReplacement <: SurveyDesign
    probs::Vector{Float64}
end

"""
Without-replacement survey design. Used to compute the Horvitz-Thompson estimate of a population total.
"""
struct WithoutReplacement <: SurveyDesign
    probs::Vector{Float64}
    joint_probs::Matrix{Float64}
end

"Simple independent survey design without replacement. In a population of size `N`, each unit is selected with equal probability."
struct SI <: SurveyDesign
    N::Int
end

"Bernoulli survey design. Individuals are included in the sample with probability `p` independently of each other."
struct Bernoulli <: SurveyDesign
    p::Float64
end

"""
Population estimate and its variance. This type is returned by the `sum` and `mean` functions
when they are applied to a `SurveyDesign`.
"""
struct SampleSum
    "Point estimate (e.g., estimated total, mean, or ratio)"
    sum::Float64
    "Variance of the estimate"
    var::Float64
end


"""
Compute a design-based estimate of a population sum using samples `xs` collected with survey design `probs`.
"""
function Base.sum(xs::AbstractVector{<:Real}, probs::SurveyDesign)
    throw(DomainError(probs, "Sample totals have not been implemented for design $(typeof(probs))."))
end

"""
Compute a design-based estimate of a nonlinear function of population totals using samples `xs` collected with survey design `probs`.
Uses Taylor series linearization with gradients from automatic differentiation.
"""
function Base.sum(f::Function, xs::Matrix{<:Real}, probs::SurveyDesign)
    throw(DomainError(probs, "Taylor estimation has not been implemented for design $(typeof(probs))."))
end

"""
Combine cluster-level estimates for multi-stage sampling designs.

This is used for cluster sampling or multi-stage designs where you have already computed
estimates within each sampled cluster. The variance properly accounts for clustering.

# Examples
```julia
# Two-stage sampling: compute within-cluster totals, then aggregate
cluster_totals = [sum(cluster1_data, SI(n1)), sum(cluster2_data, SI(n2))]
overall_total = sum(cluster_totals, SI(N_clusters))
```

# See also
R's `survey::svydesign()` with `id=~cluster`
"""
function Base.sum(xs::AbstractVector{SampleSum}, probs::SurveyDesign)
    ss = sum([x.sum for x in xs], probs)
    SampleSum(ss.sum, ss.var + cluster_variance([x.var for x in xs], probs))
end

cluster_variance(xs, probs::SI) = probs.N * mean(xs)
cluster_variance(xs, probs::WithoutReplacement) = sum(xs ./ probs.probs)


Base.:+(a::SampleSum, b::SampleSum) = SampleSum(a.sum + b.sum, a.var + b.var)
Base.:+(a::Real, b::SampleSum) = SampleSum(a + b.sum, b.var)
Base.:+(a::SampleSum, b::Real) = SampleSum(a.sum + b, a.var)

function StatsAPI.confint(a::SampleSum; level=0.95)
    α = (1 - level) / 2
    quantile.(Normal(a.sum, sqrt(a.var)), [α, 1 - α])
end

"""
Multivariate population estimates.
These structs are created when `sum` is passed Matrix and `SurveyDesign` arguments.
`SampleSums` structs are usually collapsed into a
`SampleSum` using `sum(f::Function, xs)` for stratified studies or cluster studies.

# Example
```julia
sums = sum([api_stu; enroll], SI(N))
total = sum(a -> a[1] / a[2], [sums])
```
"""
struct SampleSums
    "Estimates of population totals for each variable"
    sums_M::Vector{Float64}
    "Sample x variable matrix of observations"
    samples_nM::Matrix{Float64}
    "Population size"
    N::Int
end

Base.sum(xs::AbstractVector{<:Real}, probs::Bernoulli) =
    SampleSum(1 / probs.p * sum(xs), (1 / probs.p - 1) * sum(xs .^ 2))

Base.sum(xs::AbstractMatrix{<:Real}, probs::SI) =
    SampleSums(probs.N * vec(mean(xs; dims=1)), xs, probs.N)

function Base.sum(f::Function, xs::Vector{SampleSums})
    x0_M = sum(x.sums_M for x in xs)
    result = DiffResults.GradientResult(x0_M)
    ForwardDiff.gradient!(result, f, x0_M)
    ∇f_M = DiffResults.gradient(result)
    result_var = sum(xs; init=0.0) do x
        sum(x.samples_nM * ∇f_M, SI(x.N)).var
    end
    SampleSum(DiffResults.value(result), result_var)
end

function Base.sum(f::Function, xs::Vector{SampleSums}, probs::SI)
    N = probs.N
    x0_M = N * mean(x.sums_M for x in xs)
    result = DiffResults.GradientResult(x0_M)
    ForwardDiff.gradient!(result, f, x0_M)
    ∇f_M = DiffResults.gradient(result)
    us = [sum(x.samples_nM * ∇f_M, SI(x.N)) for x in xs]
    v = var(u.sum for u in us; corrected=true)
    SampleSum(DiffResults.value(result),
        N^2 * (1 / length(us) - 1 / N) * v + N * mean(u.var for u in us))
end

function Base.sum(xs::AbstractVector{<:Real}, probs::SI)
    N = probs.N
    if length(xs) == N
        SampleSum(sum(xs), 0.0)
    else
        m, v = mean_and_var(xs; corrected=true)
        SampleSum(N * m, N^2 * (1 / length(xs) - 1 / N) * v)
    end
end

function Base.sum(xs::AbstractVector{<:Real}, probs::WithoutReplacement)
    Δ = 1 .- (probs.probs .* probs.probs') ./ probs.joint_probs
    y = xs ./ probs.probs
    SampleSum(sum(y), y' * (Δ * y))
end

function Base.sum(f::Function, xs::Matrix{<:Real}, probs::WithoutReplacement)
    x0 = vec(sum(xs ./ reshape(probs.probs, (:, 1)); dims=1))
    result = DiffResults.GradientResult(x0)
    ForwardDiff.gradient!(result, f, x0)
    ∇f = DiffResults.gradient(result)
    u = (xs * ∇f) ./ probs.probs
    Δ = 1 .- (probs.probs .* probs.probs') ./ probs.joint_probs
    SampleSum(DiffResults.value(result), u' * (Δ * u))
end

function Base.sum(f::Function, xs::Matrix{<:Real}, probs::SI)
    N = probs.N
    x0 = N * vec(mean(xs; dims=1))
    result = DiffResults.GradientResult(x0)
    ForwardDiff.gradient!(result, f, x0)
    ∇f = DiffResults.gradient(result)
    u = xs * ∇f
    SampleSum(DiffResults.value(result), N^2 * (1 / size(xs, 1) - 1 / N) * var(u; corrected=true))
end

function GLM.lm(f::FormulaTerm, df, probs::SI)
    N = probs.N
    X, y = (modelmatrix(f, df), response(f, df))
    XX = X' * X
    β = XX \ (X'y)
    V = PDiagMat((y - X * β) .^ 2)
    n = length(y)
    SampleSum.(β, (1 - n / N) * (n / (n - 1)) * diag(XX \ (XX \ Xt_A_X(V, X))'))
end

function GLM.lm(f::FormulaTerm, df, probs::WithoutReplacement)
    X, y = (modelmatrix(f, df), response(f, df))
    Λ = Diagonal(1 ./ probs.probs)
    XX = Xt_A_X(Λ, X)
    β = XX \ (X' * Λ * y)
    R = Diagonal(y - X * β)
    Δ = 1 .- (probs.probs .* probs.probs') ./ probs.joint_probs
    V = X_A_Xt(Δ, X' * Λ * R)
    SampleSum(β, XX \ (XX \ V)')
end

pop_total(probs::SI) = probs.N
pop_total(probs::WithoutReplacement) = sum(probs.probs)

function Base.sum(f::FormulaTerm, df, df2, probs::SurveyDesign)
    X, y = (modelmatrix(f, df), response(f, df))
    X2 = modelmatrix(f, df2)
    t_x = sum(X2; dims=1)
    if f.rhs isa Tuple && f.rhs[1] isa ConstantTerm
        t_x[1] = pop_total(probs)
    end
    model_based_sum(X, t_x, y, probs)
end

function model_based_sum(X, t_x, y::AbstractVector{SampleSum}, probs::SurveyDesign)
    ss = model_based_sum(X, t_x, [yi.sum for yi in y], probs)
    SampleSum(ss.sum, ss.var + cluster_variance([yi.var for yi in y], probs))
end

# Hold on: why is the expression for WithoutReplacment so much
# simpler than for SI? What if I made g = vec(t_x * A) for SI instead?
function model_based_sum(X, t_x, y, probs::WithoutReplacement)
    Λ = Diagonal(1 ./ probs.probs)
    XX = Xt_A_X(Λ, X)
    A = XX \ (X' * Λ)
    g = vec(t_x * A)
    β = A * y
    e = y - X * β
    u = g .* e
    Δ = 1 .- (probs.probs .* probs.probs') ./ probs.joint_probs
    SampleSum(sum(g .* y), u' * (Δ * u))
end

function model_based_sum(X, t_x, y, probs::SI)
    N = probs.N
    n = length(y)
    XX = X' * X
    A = XX \ X'
    t_hat_x = N * mean(X; dims=1)
    n = size(X, 1)
    g = 1 .+ vec((n / N) * (t_x - t_hat_x) * A)
    β = A * y
    e = y - X * β
    SampleSum(N * mean(g .* y), N^2 * (1 - n / N) / n * var(g .* e; corrected=true))
end

function Base.sum(xs::AbstractVector{<:Real}, probs::WithReplacement)
    y = xs ./ probs.probs
    SampleSum(mean(y), var(y; corrected=true) / length(xs))
end

end # module Surveys
