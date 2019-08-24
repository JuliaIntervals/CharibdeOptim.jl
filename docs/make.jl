using Documenter, CharibdeOptim

makedocs(
    modules = [CharibdeOptim],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "CharibdeOptim.jl",
    authors  = "Yashvardhan Sharma",
    pages = [
        "Home" => "index.md"
         ]
)

deploydocs(
    repo = "github.com/yashcodes/CharibdeOptim.jl.git",
    target = "build",
    deps = nothing,
    make = nothing
)


