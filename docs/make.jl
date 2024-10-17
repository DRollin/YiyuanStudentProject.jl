
using Documenter, YiyuanStudentProject

DocMeta.setdocmeta!(YiyuanStudentProject, :DocTestSetup, :(using YiyuanStudentProject); recursive=true)

makedocs(;
    modules=[YiyuanStudentProject],
    authors="yiyuan.jiang@tu-braunschweig.de>, David Rollin <d.rollin@tu-braunschweig.de>",
    sitename="YiyuanStudentProject.jl",
    format=Documenter.HTML(;
        canonical="https://github.com/DRollin/YiyuanStudentProject.jl.git",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",

    ],
)

deploydocs(;
    repo="github.com/DRollin/YiyuanStudentProject.jl",
    devbranch="main",
    #push_preview=true,
)