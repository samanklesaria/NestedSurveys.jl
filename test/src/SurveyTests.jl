module SurveyTests

import Random
using NestedSurveys, RCall, DataFrames, DataFramesMeta, CSV, StatsBase, StatsModels, GLM, LinearAlgebra

function get_data()
    crime = DataFrame(CSV.File("crime_data.csv"))
    year_df = @chain crime begin
        @subset :IndexYear .<= 2011
        @subset :county .!= "STATEWIDE"
        groupby([:IndexYear, :county], sort=true)
        @combine(:rows=length(:county), :county_theft=sum(:Theft))
        @subset (:rows .> 4)
    end
    year_df[!, :total_theft] .= sum(year_df[!, :county_theft])
    innerjoin(crime, year_df, on=[:IndexYear, :county])
end

function extract_samples(g, design_sym, stage)
    N = nrow(g)
    n = min(N, 8)
    @assert n > 1
    if design_sym == :si
        ixs = sample(1:N, n; replace=false)
    elseif design_sym == :replace
        ixs = sample(1:N, n; replace=true)
    end
    prob = float(n) / N
    if stage < 2
        @transform(g[ixs, :], :fpc = N, :probs = prob, :id1 = 1:n)
    else
        @transform(g[ixs, :], :fpc2 = N, :probs2 = prob, :id2 = 1:n)
    end
end

function make_design_obj(df, design_sym, probs_col, fpc_col)
    N = df[1, fpc_col]
    if design_sym == :si
        SI(N)
    elseif design_sym == :replace
        WithReplacement(df[:, probs_col], N)
    end
end

function jl_stage0_estimate(df, design_sym, model)
    total_theft = df[1, :total_theft]
    design = make_design_obj(df, design_sym, :probs, :fpc)
    if model == :sum
        sum(df[!,:Burglary], design)
    elseif model == :regress
        sum(@formula(Burglary ~ 1 + Theft), df, (; Theft=[total_theft]), design)
    end
end

function r_stage0_estimate(df, design_sym, model)
    if design_sym == :si
        R"design <- svydesign(id=~1, fpc=~fpc, data=df)"
    elseif design_sym == :replace
        R"design <- svydesign(id=~1, probs=~probs, data=df)"
    end
    if model == :sum
        SampleSum(rcopy(R"with_var(svytotal(~Burglary, design))")...)
    elseif model == :regress
        totals = select(df, :fpc => "(Intercept)", :total_theft=>:Theft)[1:1,:]
        @rput totals
        R"design <- calibrate(design, formula=~Theft, population=totals, calfun='linear')"
        SampleSum(rcopy(R"with_var(svytotal(~Burglary, design))")...)
    end
end

function test_stage0(df, design_sym, model)
    println("d-$design_sym m-$model")
    df = extract_samples(df, design_sym, 1)
    @rput df
    r = r_stage0_estimate(df, design_sym, model)
    jl = jl_stage0_estimate(df, design_sym, model)
    @assert jl ≈ r "$jl ≉ $r"
end

function testit()
    Random.seed!(7)
    R"""
    library(survey)
    library(tidyverse)
    with_var <- function(result) {
        c(result[[1]], SE(result)^2[[1]])
    }
    """
    df = get_data()
    df = @subset(df, :IndexYear .== 2004)
    for design_sym in [:si, :replace]
        for model in [:sum, :regress]
            test_stage0(df, design_sym, model)
        end
    end
end

end
