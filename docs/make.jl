using PGA2D
using Documenter

DocMeta.setdocmeta!(PGA2D, :DocTestSetup, :(using PGA2D); recursive=true)

makedocs(;
    modules=[PGA2D],
    authors="Andreas Tell <atell@soundtheory.com> and contributors",
    repo="https://github.com/ATell-SoundTheory/PGA2D.jl/blob/{commit}{path}#{line}",
    sitename="PGA2D.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ATell-SoundTheory.github.io/PGA2D.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Triangle centers" => "triangle_centers.md",
    ],
)

deploydocs(;
    repo="github.com/ATell-SoundTheory/PGA2D.jl",
)
