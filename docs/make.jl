using Multibody
using Documenter
using CairoMakie
ENV["JULIA_DEBUG"]=Documenter # Enable this for debugging

DocMeta.setdocmeta!(Multibody, :DocTestSetup, :(using Multibody); recursive = true)

makedocs(;
         modules = [Multibody],
         authors = "JuliaHub Inc.",
         #  strict = [:example_block, :setup_block, :eval_block],
         sitename = "Multibody.jl",
         warnonly = [:missing_docs, :cross_references, :example_block, :docs_block],
         pagesonly = true,
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", nothing) == "true",
                                  edit_link = nothing),
         pages = [
             "Home" => "index.md",
             "Tutorials" => [
                 "Getting started: Pendulum" => "examples/pendulum.md",
             ],
             "Examples" => [
                 "Spring-damper system" => "examples/spring_damper_system.md",
                 "Spring-mass system" => "examples/spring_mass_system.md",
                 "Three springs (series forces)" => "examples/three_springs.md",
                 "Sensors" => "examples/sensors.md",
                 "Spherical pendulum" => "examples/spherical_pendulum.md",
                 "Gearbox" => "examples/gearbox.md",
                 "Free motions" => "examples/free_motion.md",
                 "Kinematic loops" => "examples/kinematic_loops.md",
                 "Industrial robot" => "examples/robot.md",
             ],
             "3D rendering" => "rendering.md",
         ])

deploydocs(;
           repo = "github.com/JuliaComputing/Multibody.jl",
           devbranch = "main")
