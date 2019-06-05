using Documenter, CharibdeOptim

makedocs(;
    modules=[CharibdeOptim],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/yashcodes/CharibdeOptim.jl/blob/{commit}{path}#L{line}",
    sitename="CharibdeOptim.jl",
    authors="Chris de Graaf, Invenia Technical Computing Corporation",
    assets=String[],
)

deploydocs(;
    repo="github.com/yashcodes/CharibdeOptim.jl",
)
