# NestedSurveys.jl

A Julia package for design-based inference in survey sampling. This package provides functionality comparable to R's `survey` package, with native Julia performance and integration with DataFrames.jl. It is similar in spirit to the [Survey.jl](https://xkdr.github.io/Survey.jl) package, but allows for true multi-stage sampling with different designs at each stage.

## Overview

`NestedSurveys.jl` implements methods for analyzing data from complex survey designs, including:

- Simple random sampling (SRS)
- Stratified sampling
- One-stage cluster sampling
- Two-stage cluster sampling
- Taylor series variance estimation
- Confidence intervals
- Ratio estimation

The main exported type is `SampleSum`, which stores both an estimate and its variance. Sampling designs are specified with `SurveyDesign` structs, including `SI` (simple random sampling without replacement), `WithReplacement`,  `WithoutReplacement`, and `Bernoulli`.
    
## Simple Random Sampling

For simple random sampling without replacement, use `sum` with an `SI` object that specifies the population size.

**Julia:**
```julia
using NestedSurveys, DataFramesMeta
result = @combine(apisrs, :total = sum(:enroll, SI(N)))
```

**R equivalent:**
```r
library(survey)
srs_design <- svydesign(id=~1, fpc=~fpc, data=apisrs)
svytotal(~enroll, srs_design)
```

The `sum` function computes the Horvitz-Thompson estimator for the total and its variance, accounting for the finite population correction.

## Stratified Sampling

For stratified sampling, compute subtotals within each stratum, then combine them.

**Julia:**
```julia
strat_result = @chain apistrat begin
    @groupby(:stype)
    @combine(:subtotal = sum(:enroll, SI(Int(:fpc[1]))))
    @combine(:total = sum(:subtotal))
end
```

**R equivalent:**
```r
strat_design <- svydesign(id=~1, fpc=~fpc, strata=~stype, data=apistrat)
svytotal(~enroll, strat_design)
```

The `SampleSum` type supports addition, so stratified estimates can be combined by summing the `subtotal` column.

## One-Stage Cluster Sampling

For one-stage cluster sampling, first aggregate within clusters, then use `sum` with `SI` on the cluster totals.

**Julia:**
```julia
gdf = groupby(cal_crime, :county)
@chain gdf begin
    @combine(:subtotal = sum(:Burglary))
    @combine(:total = sum(:subtotal, SI(N_counties)))
end
```

**R equivalent:**
```r
stage1_design <- svydesign(id=~county, fpc=~fpc, data=cal_crime)
svytotal(~Burglary, stage1_design)
```

The key is to first compute totals within each sampled cluster, then treat those cluster totals as the observations for `sum`.

## Two-Stage Cluster Sampling

For two-stage sampling, apply `sum` with `SI` twice: once within primary sampling units (PSUs), then across PSUs.

**Julia:**
```julia
@chain df begin
    @groupby(:county)
    @combine(:subtotal = sum(:Burglary, SI(county_sizes[first(:county)])))
    @combine(:total = sum(:subtotal, SI(N_counties)))
end
```

**R equivalent:**
```r
stage2_design <- svydesign(id=~county+id, fpc=~fpc+fpc2, data=df)
svytotal(~Burglary, stage2_design)
```

The variance calculation properly accounts for both stages of sampling through the nested application of `sum`.

## Ratio Estimation

For ratio estimation or other nonlinear functions of totals, use `taylor` with a `Function` argument for Taylor series linearization.

**Julia:**
```julia
# Estimate ratio of api.stu to enroll
ratio_result = @combine(apisrs, :total = 
    taylor(a -> a[1] / a[2], g-> sum(g([:api_stu :enroll]), SI(Int(:fpc[1])))))
```

**R equivalent:**
```r
srs_design <- svydesign(id=~1, fpc=~fpc, data=apisrs)
svyratio(~api.stu, ~enroll, srs_design)
```

The `sum` function uses automatic differentiation to compute the Taylor series approximation to the variance of nonlinear estimators.

## Regression-Assisted Estimation

For regression-assisted (calibration) estimation, use `sum` with a formula, sample data, population data, and design.

**Julia:**
```julia
assisted_result = sum(@formula(api_stu ~ 1 + enroll), apisrs, (; enroll=[4e6]), SI(6194))
```

**R equivalent:**
```r
srs_design <- svydesign(id=~1, fpc=~fpc, data=apisrs)
svytotal(~api.stu, calibrate(srs_design, ~enroll, c('(Intercept)'=6194, enroll=4e6)))
```

## Sampling With Replacement

For sampling with replacement, use `sum` with a `WithReplacement` object containing sampling probabilities.

```julia
result = sum(observations, WithReplacement(sample_probs))
```

This computes the Hansen-Hurwitz estimator for the total.

## Unequal Probability Sampling Without Replacement

For unequal probability sampling without replacement, use `sum` with a `WithoutReplacement` object containing inclusion probabilities and joint inclusion probabilities.

```julia
result = sum(observations, WithoutReplacement(inclusion_probs, joint_inclusion_probs))
```

This computes the Horvitz-Thompson estimator with arbitrary inclusion probabilities.

## API Reference

```@autodocs
Modules = [NestedSurveys]
```
