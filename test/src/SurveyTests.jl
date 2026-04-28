module SurveyTests

import Random
using NestedSurveys, RCall, DataFrames, DataFramesMeta, CSV, StatsBase, StatsModels, GLM, LinearAlgebra

function simple_ratio_test()
    X = randn(5, 2)
    taylor(y -> y[1] / y[2]) do g
        sum(g(X), SI(50))
    end
end


function get_data()
    crime = DataFrame(CSV.File("crime_data.csv"))
    year_df = @chain crime begin
        @subset :IndexYear .<= 2011
        @subset :county .!= "STATEWIDE"
        groupby([:IndexYear, :county], sort=true)
        combine(nrow => :rows)
        @subset (:rows .> 4)
    end
    innerjoin(crime, year_df, on=[:IndexYear, :county])
end

function testit()
    Random.seed!(7)
    R"""
    library(survey)
    library(tidyverse)
    with_var <- function(result) {
        c(result[[1]], SE(result)^2[[1]])
    }
    pop_totals <- function(df) {
      df |>
        pivot_longer(c(fpc, Theft), names_to = "col") |>
        mutate(name = if_else(col == "fpc",
          paste0("factor(IndexYear)", IndexYear),
          paste0("factor(IndexYear)", IndexYear, ":Theft"))) |>
        select(name, value) |>
        deframe()
    }
    """
    # @testset "cal_crime" begin
    crime = get_data()
    totals = @by crime :IndexYear :county begin
        :Theft = sum(:Theft)
    end
    for clusters in [0, 1, 2]
        for stratified in [false, true]
            for design in [:si, :pps, :replace]
                for model1 in [:sum, :regress, :ratio] # mean
                    if clusters == 2
                        design2_options = [:si, :pps, :replace]
                        if model1 == :ratio
                            model2_options = [:sum]
                        else
                            model2_options = [:sum, :regress]
                        end
                    else
                        design2_options = [:si]
                        model2_options = [:sum]
                    end
                    for design2 in design2_options
                        for model2 in model2_options
                            println("c-$clusters s-$stratified d1-$design d2-$design2 m1-$model1 m2-$model2")
                            df, t = make_df(crime, totals, clusters, stratified, design, design2)
                            jl = jl_estimate(df, stratified, clusters, model1, model2, t, design, design2)
                            r = r_estimate(df, stratified, clusters, model1, model2, t, design, design2)
                            if !isnothing(r)
                                @assert jl ≈ r "$jl ≉ $r"
                            end
                        end
                    end
                end
            end
        end
    end
end

function extract_samples(g, design_sym, stage)
    N = nrow(g)
    n = min(nrow(g), 8)
    @assert n > 1
    if design_sym == :si
        ixs = sample(1:N, n; replace=false)
        prob = float(n) / N
    elseif design_sym == :pps
        ixs = sample(1:N, n; replace=false)
        prob = float(n) / N
    elseif design_sym == :replace
        ixs = sample(1:N, n; replace=true)
        prob = 1.0 / N
    end
    if stage < 2
        @transform(g[ixs, :], :fpc = N, :probs = prob, :id1 = 1:n)
    else
        @transform(g[ixs, :], :fpc2 = N, :probs2 = prob, :id2 = 1:n)
    end
end

"""
Returns both sample data and population totals.
Under stratification, the totals are grouped by :IndexYear
Sample data has fields:
- :id1 (PSU)
- :IndexYear
- :id2 (SSU if stages == 2)
- :fpc (number of PSUs)
- :fpc2 (number of SSUs)
- :probs (PSU inclusion probs)
- :probs2 (SSU inclusion probs)
- :Theft
- :Burglary
"""
function make_df(crime, totals, stages, stratified, design, design2)
    function apply_clusters(df)
        if stages == 0
            return extract_samples(df, design, 1)
        end
        counties = @chain df begin
            groupby(:county)
            @combine(:Theft = sum(:Theft))
            extract_samples(design, 1)
            @select(:county, :fpc, :probs, :id1)
        end
        df = innerjoin(counties, df, on=:county)
        if stages == 1
            df
        else
            @chain df begin
                @groupby(:id1)
                combine(x -> extract_samples(x, design2, 2))
            end
        end
    end

    if stratified
        combine(groupby(crime, :IndexYear), apply_clusters), groupby(totals, :IndexYear)
    else
        apply_clusters(@subset(crime, :IndexYear .== 2004)), @subset(totals, :IndexYear .== 2004)
    end
end

function r_estimate(df, stratified, clusters, model1, model2, totals, design, design2)
    if model2 != :sum
        return nothing # TODO
    end
    @rput df
    id_arg = clusters == 0 ? "id=~1" : clusters == 1 ? "id=~id1" : "id=~id1+id2"
    if design == :si && design2 == :si
        fpc_arg = clusters == 2 ? "fpc=~fpc+fpc2" : "fpc=~fpc"
    elseif design == :replace && design2 == :replace
        fpc_arg = clusters == 2 ? "probs=~probs+probs2" : "probs=~probs"
    else
        return nothing # Not supported in R
    end
    strata_arg = stratified ? ", strata=~IndexYear, nest=TRUE" : ""
    design_str = "design_r <- svydesign($id_arg, $fpc_arg$strata_arg, data=df)"
    reval(design_str)
    if model1 == :sum
        SampleSum(rcopy(R"with_var(svytotal(~Burglary, design_r))")...)
    elseif model1 == :ratio
        if clusters > 1
            return nothing # TODO
        end
        SampleSum(rcopy(R"with_var(svyratio(~Burglary, ~Theft, design_r))")...)
    elseif model1 == :regress
        r_totals = @combine(totals, :Theft = sum(:Theft))
        if stratified
            counts = @by(df, :IndexYear, :fpc = first(:fpc))
            r_totals = innerjoin(r_totals, counts, on=:IndexYear)
        else
            r_totals = insertcols(r_totals, "(Intercept)" => df[1, :fpc])
        end

        @rput r_totals
        if stratified
            R"r_totals <- pop_totals(r_totals)"
            formula="formula=~0+factor(IndexYear)+factor(IndexYear):Theft"
        else
            formula = "formula=~Theft"
        end
        if clusters != 0
           return nothing # TODO
        end
        reval("c_design <- calibrate(design_r, $formula, population=r_totals, calfun='linear')")
        SampleSum(rcopy(R"with_var(svytotal(~Burglary, c_design))")...)
    end
end

function make_design_obj(df, design_sym, probs_col, fpc_col)
    N = df[1, fpc_col]
    if design_sym == :si
        SI(N)
    elseif design_sym == :replace
        WithReplacement(df[:, probs_col], N)
    elseif design_sym == :pps
        n = nrow(df)
        probs = df[:, probs_col]
        joint = n * (n - 1) / (N * (N - 1))
        joint_mat = fill(joint, (n, n))
        joint_mat[diag(CartesianIndices(joint_mat))] .= probs
        WithoutReplacement(probs, joint_mat)
    end
end

function cluster_totals(df, design_sym, probs_col, fpc_col, model, totals, g)
    d = make_design_obj(df, design_sym, probs_col, fpc_col)
    if g == nothing
        total_theft = sum(totals[!, :Theft])
        if model == :sum
            @combine(df,
                :fpc = first(:fpc),
                :probs = first(:probs), # for if :probs_col = :probs2
                :Theft = total_theft,
                :Burglary = sum(:Burglary, d))
        elseif model == :regress
            total = sum(@formula(Burglary ~ 1 + Theft), df, totals, d)
            DataFrame(:fpc => df[1, :fpc], :probs => df[1, :probs],
                :Theft => total_theft, :Burglary => total)
        end
    else
        @combine(df, :Burglary = Ref(sum(g([:Burglary :Theft]), d)))
    end
end

function jl_estimate(df, stratified, clusters, model1, model2, totals, design, design2)
    function get_subtotals(df, totals, g)
        if clusters == 0
            cluster_totals(df, design, :probs, :fpc, model1, totals, g)
        elseif clusters == 1
            @chain df begin
                groupby(:id1)
                @combine(:Burglary = sum(:Burglary),
                    :Theft = sum(:Theft),
                    :fpc = first(:fpc), :probs = first(:probs))
                cluster_totals(design, :probs, :fpc, model1, totals, g)
            end
        else
            g_totals = groupby(totals, :county)
            @chain df begin
                groupby(:id1)
                combine(x -> cluster_totals(x, design2, :probs2, :fpc2,
                        model2, g_totals[(; county=x[1, :county])], nothing))
                cluster_totals(design, :probs, :fpc, model1, totals, g)
            end
        end
    end

    function calculate(g)
        if stratified
            tot = @chain df begin
                groupby(:IndexYear)
                combine(df ->
                    get_subtotals(df, totals[(; IndexYear=df[1, :IndexYear])], g))
                @combine(:Burglary = Ref(sum(:Burglary)))
            end
        else
            tot = get_subtotals(df, totals, g)
        end
        tot[1, :Burglary]
    end

    if model1 == :ratio
        taylor(calculate, y -> y[1] / y[2])
    else
        calculate(nothing)
    end
end

end
