# Surveys.jl

A Julia package for design-based inference in survey sampling. This package provides functionality comparable to R's `survey` package, with native Julia performance and integration with DataFrames.jl. [![][docs-dev-img]][docs-dev-url]

## Overview

`Surveys.jl` implements methods for analyzing data from complex survey designs, including:

- Simple random sampling (SRS)
- Stratified sampling
- One-stage cluster sampling
- Two-stage cluster sampling
- Taylor series variance estimation
- Ratio estimation
- Regression-assisted estimation

The main exported type is `SampleSum`, which stores both an estimate and its variance. Sampling designs are specified using `SI` (simple random sampling without replacement), `WithReplacement`, and `WithoutReplacement` structs.

## Simple Random Sampling

For simple random sampling without replacement, use `sum` with an `SI` object that specifies the population size.

```julia
using Surveys, DataFramesMeta
result = @combine(apisrs, :total = sum(:enroll, SI(N)))
```

The `sum` function computes the Horvitz-Thompson estimator for the total and its variance, accounting for the finite population correction.

## Stratified Sampling

For stratified sampling, compute subtotals within each stratum, then combine them.

```julia
strat_result = @chain apistrat begin
    @groupby(:stype)
    @combine(:subtotal = sum(:enroll, SI(Int(:fpc[1]))))
    @combine(:total = sum(:subtotal))
end
```

The `SampleSum` type supports addition, so stratified estimates can be combined by summing the `subtotal` column.

## One-Stage Cluster Sampling

For one-stage cluster sampling, first aggregate within clusters, then use `sum` with `SI` on the cluster totals.

```julia
gdf = groupby(cal_crime, :county)
@chain gdf begin
    @combine(:subtotal = sum(:Burglary))
    @combine(:total = sum(:subtotal, SI(N_counties)))
end
```

## Two-Stage Cluster Sampling

For two-stage sampling, apply `sum` with `SI` twice: once within primary sampling units (PSUs), then across PSUs.

```julia
@chain df begin
    @groupby(:county)
    @combine(:subtotal = sum(:Burglary, SI(county_sizes[first(:county)])))
    @combine(:total = sum(:subtotal, SI(N_counties)))
end
```

The variance calculation properly accounts for both stages of sampling through the nested application of `sum`.

## Ratio Estimation

For ratio estimation or other nonlinear functions of totals, use `sum` with a `Function` argument for Taylor series linearization.

```julia
# Estimate ratio of api.stu to enroll
ratio_result = @combine(apisrs, :total = 
    sum((a -> a[1] / a[2]), [:api_stu :enroll], SI(Int(:fpc[1]))))
```

The `sum` function uses automatic differentiation to compute the Taylor series approximation to the variance of nonlinear estimators.

## Linearization with Stratification and Clustering

When `sum` is passed a Matrix and `SI` object, it creates a `SampleSums` object. This can be passed to `sum` to get clustered or stratified Taylor series variance estimates.

```julia
@chain apistrat begin
    @groupby(:stype)
    @combine(:subtotal = sum([:api_stu :enroll], SI(Int(:fpc[1]))))
    @combine(:total = sum(a->a[1] / a[2], :subtotal))
end
```

## Coefficient Estimation

For regression coefficient estimation with design-based variance, use `π_lm` with a formula and design specification.

```julia
π_lm(@formula(api_stu ~ 1 + enroll), apisrs, SI(Int(apisrs[1, :fpc])))
```

The `π_lm` function returns a vector of `SampleSum` objects, one for each coefficient, with design-based variance estimates.

## Regression-Assisted Estimation

For regression-assisted (calibration) estimation, use `sum` with a formula, sample data, population data, and design.

```julia
assisted_result = sum(@formula(api_stu ~ 1 + enroll), apisrs, (; enroll=[4e6]), SI(Int(apisrs[1, :fpc])))
```

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://samanklesaria.github.io/Surveys.jl