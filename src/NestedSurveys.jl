module NestedSurveys
using Statistics, StatsBase, StatsAPI, StatsModels, DiffResults, ForwardDiff, PDMats
using LinearAlgebra, Distributions
using SizeCheck
include("docboilerplate.jl")

export SampleSum, SI, WithReplacement, WithoutReplacement, SurveyDesign, taylor

# Abstract Interfaces

"""
Abstract Survey Design. When passed as an argument to `sum` or `mean`, computes the design-based estimate and variance.
"""
abstract type SurveyDesign end

"""
With-replacement survey design. Used for probability-weighted ratio estimators (Hansen-Hurwitz) for population totals and means.
"""
Base.@kwdef struct WithReplacement <: SurveyDesign
    "Marginal sample inclusion probabilities"
    probs::Vector{Float64}
    "Population size"
    N::Union{Nothing,Int} = nothing
end

"""
Without-replacement survey design. Used to compute the Horvitz-Thompson estimate of a population total.
"""
struct WithoutReplacement <: SurveyDesign
    "Marginal sample inclusion probabilities"
    probs::Vector{Float64}
    "Pairwise joint sample inclusion probabilities"
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

cluster_estimate(x::Real) = x
cluster_estimate(x::SampleSum) = x.sum

cluster_var(x::Real) = 0.0
cluster_var(x::SampleSum) = x.var

Base.:≈(a::SampleSum, b::SampleSum) = (a.sum ≈ b.sum) && (a.var ≈ b.var)

Base.:/(a::SampleSum, b::Real) = SampleSum(a.sum / b, a.var / b^2)
Base.:*(a::SampleSum, b::Real) = SampleSum(a.sum * b, a.var * b^2)
Base.:*(a::Real, b::SampleSum) = SampleSum(b.sum * a, b.var * a^2)

"""
Combine cluster-level estimates for multi-stage sampling designs. This is used for cluster sampling or
multi-stage designs where you have already computed estimates within each sampled cluster.
The returned variance properly accounts for clustering.
"""
function Base.sum(xs::AbstractVector, design::SurveyDesign)
    ss = sum(cluster_estimate.(xs), design)
    SampleSum(ss.sum, ss.var + sum(cluster_var.(xs), design).sum)
end

function Base.sum(xs::AbstractMatrix, design::SurveyDesign)
    mapslices(xs, dims=1) do x
        sum(x, design)
    end
end

"""
Compute a design-based estimate of a population mean using samples `xs` collected with survey design `probs`.
"""
Statistics.mean(xs::AbstractVector, probs::SurveyDesign) =
    taylor(y -> y[1] / y[2]) do g
        sum(g([xs ones(length(xs))]), probs)
    end

function taylor(sum_of::Function, f::Function)
    xs = vec(cluster_estimate.(sum_of(x -> cluster_estimate.(x))))
    result = DiffResults.GradientResult(xs)
    ForwardDiff.gradient!(result, f, xs)
    ∇f_M = DiffResults.gradient(result)
    var_of_E = sum_of(x -> cluster_estimate.(x) * ∇f_M).var
    E_of_var = sum_of(x -> cluster_var.(x) * ∇f_M).sum
    SampleSum(DiffResults.value(result), var_of_E + E_of_var)
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
    SampleSum(ss.sum, ss.var + sum([x.var for x in xs], probs).sum)
end

Base.:+(a::SampleSum, b::SampleSum) = SampleSum(a.sum + b.sum, a.var + b.var)
Base.:+(a::Real, b::SampleSum) = SampleSum(a + b.sum, b.var)
Base.:+(a::SampleSum, b::Real) = SampleSum(a.sum + b, a.var)

function StatsAPI.confint(a::SampleSum; level=0.95)
    α = (1 - level) / 2
    quantile.(Normal(a.sum, sqrt(a.var)), [α, 1 - α])
end

Base.sum(xs::AbstractVector{<:Real}, probs::Bernoulli) =
    SampleSum(1 / probs.p * sum(xs), (1 / probs.p - 1) * sum(xs .^ 2))

function Base.sum(xs::AbstractVector{<:Real}, probs::WithReplacement)
    y = xs ./ probs.probs
    SampleSum(mean(y), var(y; corrected=true) / length(xs))
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

pop_total(probs::SI) = probs.N
pop_total(probs::WithoutReplacement) = sum(probs.probs)

function pop_total(probs::WithReplacement)
    if isnothing(probs.N)
        error("Population total unknown; cannot use intercept in model based estimate")
    end
    probs.N
end

function Base.sum(f::FormulaTerm, df, df2, probs::SurveyDesign)
    y = modelcols(f.lhs, df)
    X = modelmatrix(f.rhs, df)
    X2 = modelmatrix(f.rhs, df2)
    s = sum(X2; dims=1)
    if f.rhs isa Tuple && f.rhs[1] isa ConstantTerm
        s[1] = pop_total(probs)
    end
    model_based_sum(X, s, y, probs)
end

function model_based_sum(X, s, y::AbstractVector{SampleSum}, probs::SurveyDesign)
    ss = model_based_sum(X, s, cluster_estimate.(y), probs)
    SampleSum(ss.sum, ss.var + sum([yi.var for yi in y], probs).sum)
end

function T_estimator(X, probs::T) where {T<:Union{WithoutReplacement,WithReplacement}}
    Xt_A_X(PDiagMat(1 ./ probs.probs), X)
end

function T_estimator(X, probs::SI)
    n = size(X, 1)
    Xt_A_X(ScalMat(n, probs.N / n), X)
end

@sizecheck function t_estimator(X_MD, y_M, probs::T) where {T<:Union{WithoutReplacement,WithReplacement}}
    X_MD' * (y_M ./ probs.probs)
end

@sizecheck function t_estimator(X_MD, y_M, probs::SI)
    (probs.N ./ M) * (X_MD' * y_M)
end

@sizecheck function model_based_sum(X_ND, s_1D, y_N, probs::SurveyDesign)
    @assert N > 1
    T_DD = T_estimator(X_ND, probs)
    t_D = t_estimator(X_ND, y_N, probs)
    β_D = T_DD \ t_D
    g_N = vec(s_1D * (T_DD \ X_ND')) .* (y_N - X_ND * β_D)
    SampleSum(dot(vec(s_1D), β_D), sum(g_N, probs).var)
end

end # module Surveys
