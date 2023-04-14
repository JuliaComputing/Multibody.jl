using Multibody
using Documenter

DocMeta.setdocmeta!(Multibody, :DocTestSetup, :(using Multibody); recursive = true)

makedocs(;
         modules = [Multibody],
         authors = "Yingbo Ma <mayingbo5@gmail.com> and contributors",
         sitename = "Multibody.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", nothing) == "true",
                                  edit_link = nothing),
         pages = [
             "Home" => "index.md",
         ])

deploydocs(;
           repo = "github.com/JuliaComputing/Multibody.jl",
           devbranch = "master")
