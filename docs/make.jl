using Documenter

repo_url = "github.com/control-toolbox/GeometricPreconditioner.jl"

makedocs(;
    warnonly=:cross_references,
    sitename="Geometric preconditioner",
    format=Documenter.HTML(;
        repolink="https://" * repo_url,
        prettyurls=false,
        size_threshold_ignore=["2D-example.md"],
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
    ),
    pages=[
        "Introduction" => "index.md",
        "Indirect shooting" => "2D-example.md",
        "Geometric preconditioner" => "2D-preconditioner.md",
    ],
)

deploydocs(; repo=repo_url * ".git", devbranch="main")
