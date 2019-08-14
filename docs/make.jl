using Documenter, CharibdeOptim

makedocs(
    modules = [CharibdeOptim],
    doctest = true,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages=["Home" => "index.md"],
    sitename = "CharibdeOptim.jl",
    authors = "Yashvardhan Sharma",
)

deploydocs(;
    repo="github.com/yashcodes/CharibdeOptim.jl",
    target = "build",
    deps = nothing,
    make = nothing,
)
