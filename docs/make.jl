using Multibody
using Documenter

DocMeta.setdocmeta!(Multibody, :DocTestSetup, :(using Multibody); recursive = true)

makedocs(;
         modules = [Multibody],
         authors = "Yingbo Ma <mayingbo5@gmail.com> and contributors",
         repo = "https://github.com/YingboMa/Multibody.jl/blob/{commit}{path}#{line}",
         sitename = "Multibody.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://YingboMa.github.io/Multibody.jl",
                                  edit_link = "master",
                                  assets = String[]),
         pages = [
             "Home" => "index.md",
             "Tutorials" => [
                 "Hello world: Pendulum" => "examples/pendulum.md",
             ],
         ])

deploydocs(;
           repo = "github.com/JuliaComputing/Multibody.jl",
           devbranch = "master")
