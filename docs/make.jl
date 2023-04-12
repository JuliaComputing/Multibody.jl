using Multibody
using Documenter

DocMeta.setdocmeta!(Multibody, :DocTestSetup, :(using Multibody); recursive = true)

makedocs(;
         modules = [Multibody],
         authors = "Yingbo Ma <mayingbo5@gmail.com> and contributors",
         repo = "https://github.com/YingboMa/Multibody.jl/blob/{commit}{path}#{line}",
         sitename = "Multibody.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", nothing) == "true",
                                  edit_link = nothing),
         pages = [
             "Home" => "index.md",
         ])

deploydocs(;
           repo = "github.com/YingboMa/Multibody.jl",
           devbranch = "master")
