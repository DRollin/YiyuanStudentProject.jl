
using Documenter, YiyuanStudentProject

DocMeta.setdocmeta!(YiyuanStudentProject, :DocTestSetup, :(using YiyuanStudentProject); recursive=true)

include("generate.jl")

makedocs(;
    modules=[YiyuanStudentProject],
    authors="yiyuan.jiang@tu-braunschweig.de>, David Rollin <d.rollin@tu-braunschweig.de>",
    sitename="Chemo-mechanical.jl",
    format=Documenter.HTML(;
        canonical="https://github.com/DRollin/YiyuanStudentProject.jl.git",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Documentation" => [
            "Documentation overview" => "documentation/index.md",
            "documentation/fine_scale.md",
            "documentation/upscaling.md",
            "documentation/macro_scale.md",
        ],

        "Reference" => [
            "Reference overview" => "reference/index.md",
            "reference/fine_scale.md",
            "reference/upscaling.md",
            "reference/macro_scale.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/DRollin/YiyuanStudentProject.jl",
    devbranch="main",
    push_preview=true,
)