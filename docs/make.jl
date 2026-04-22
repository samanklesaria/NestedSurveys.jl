using Documenter
using NestedSurveys

makedocs(
    sitename="NestedSurveys",
    format=Documenter.HTML(sidebar_sitename=false),
    pages=[],
    modules=[NestedSurveys]
)

deploydocs(
    repo="github.com/samanklesaria/NestedSurveys.jl.git",
    devbranch="main"
)
