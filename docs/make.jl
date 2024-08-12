using Documenter

makedocs(
    warnonly = :cross_references,
    sitename = "Geometric preconditioner",
    format = Documenter.HTML(
        prettyurls = false, 
        #size_threshold_ignore = [""],
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
    ),
    pages = [
        "Introduction" => "index.md",
        "2D example"   => "2D-example.md",
    ]
)

deploydocs(
    repo = "github.com/control-toolbox/preconditioning.git",
    devbranch = "main"
)