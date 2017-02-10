using Documenter, CutPruners

makedocs()

deploydocs(
    deps = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/JuliaPolyhedra/CutPruners.jl.git",
    julia  = "release"
)
