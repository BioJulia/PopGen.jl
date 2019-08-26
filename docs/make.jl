using Documenter, PopGen

makedocs(;
    modules=[PopGen],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/pdimens/PopGen.jl/blob/{commit}{path}#L{line}",
    sitename="PopGen",
    authors="Pavel Dimens, Jason Selwyn",
    assets=String[],
)

deploydocs(;
    repo="github.com/pdimens/PopGen.jl",
)
