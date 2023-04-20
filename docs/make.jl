using Multibody
using Documenter

DocMeta.setdocmeta!(Multibody, :DocTestSetup, :(using Multibody); recursive = true)

makedocs(;
         modules = [Multibody],
         authors = "JuliaHub Inc.",
        #  strict = [:example_block, :setup_block, :eval_block],
         sitename = "Multibody.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", nothing) == "true",
                                  edit_link = nothing),
         pages = [
             "Home" => "index.md",
             "Tutorials" => [
                 "Hello world: Pendulum" => "examples/pendulum.md",
             ],
         ])

deploydocs(;
           repo = "github.com/JuliaComputing/Multibody.jl",
           devbranch = "main")
