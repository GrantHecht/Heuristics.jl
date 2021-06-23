using Heuristics
using Documenter

DocMeta.setdocmeta!(Heuristics, :DocTestSetup, :(using Heuristics); recursive=true)

makedocs(;
    modules=[Heuristics],
    authors="Grant Hecht",
    repo="https://github.com/GrantHecht/Heuristics.jl/blob/{commit}{path}#{line}",
    sitename="Heuristics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://GrantHecht.github.io/Heuristics.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/GrantHecht/Heuristics.jl",
    devbranch="main",
)
