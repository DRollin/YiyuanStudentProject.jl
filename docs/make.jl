using TimerOutputs

dto = TimerOutput()
reset_timer!(dto)

using Documenter, YiyuanStudentProject, DocumenterCitations

DocMeta.setdocmeta!(YiyuanStudentProject, :DocTestSetup, :(using YiyuanStudentProject); recursive=true)

include("generate.jl")

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "assets", "refs.bib"),
    style=:numeric  # default
)

makedocs(;
    modules=[YiyuanStudentProject],
    authors="yiyuan.jiang@tu-braunschweig.de>, David Rollin <d.rollin@tu-braunschweig.de>",
    sitename="Chemo-mechanical.jl",
    format=Documenter.HTML(;
        canonical="https://github.com/DRollin/YiyuanStudentProject.jl.git",
        edit_link="main",
        assets=String["assets/citations.css"],
    ),
    pages=[
        "Home" => "index.md",
        "Documentation" => [
            "Documentation overview" => "documentation/index.md",
            "documentation/fine_scale.md",
            "documentation/multi_scale.md",
            "documentation/reference.md",
        ],
        "Example"  => "examples/example_1.md",
        "Reference" => [
            "Reference overview" => "reference/index.md",
            "reference/types.md",
            "reference/sub_scale.md",
            "reference/upscaling.md",
            "reference/macro_scale.md",
        ],
    ],
    plugins=[bib],
)

deploydocs(;
    repo="github.com/DRollin/YiyuanStudentProject.jl",
    devbranch="main",
    push_preview=true,
)