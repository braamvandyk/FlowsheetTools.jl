using FlowsheetTools
using Documenter

DocMeta.setdocmeta!(FlowsheetTools, :DocTestSetup, :(using FlowsheetTools); recursive=true)

makedocs(;
    modules=[FlowsheetTools],
    authors="Braam van Dyk <braam.vandyk@gmail.com> and contributors",
    repo="https://github.com/braamvandyk/FlowsheetTools.jl/blob/{commit}{path}#{line}",
    sitename="FlowsheetTools.jl",
    format=Documenter.HTML(;
        # prettyurls=get(ENV, "CI", "false") == "true",
        prettyurls=true,
        canonical="https://braamvandyk.github.io/FlowsheetTools.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/UserName/FlowsheetTools.jl",
    devbranch="main",
)
