using Multibody
using Documenter
using GLMakie
GLMakie.activate!()

ENV["JULIA_DEBUG"]=Documenter # Enable this for debugging
ENV["DOCS_BUILD"] = true # used to lower the default frame rate in animations for the docs

DocMeta.setdocmeta!(Multibody, :DocTestSetup, :(using Multibody); recursive = true)

makedocs(;
         modules = [Multibody],
         authors = "JuliaHub Inc.",
         #  strict = [:example_block, :setup_block, :eval_block],
        #  remotes = Dict(
        #     dirname(dirname(pathof(Multibody))) => (Remotes.GitHub("JuliaComputing", "Multibody.jl"), "0"),
        # ),
         sitename = "Multibody.jl",
         warnonly = [:missing_docs, :cross_references, :docs_block],
         pagesonly = true,
        #  draft = true,
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
                "Prescribed motions" => "examples/prescribed_pose.md",
                "Kinematic loops" => "examples/kinematic_loops.md",
                "Industrial robot" => "examples/robot.md",
                "Ropes, cables and chains" => "examples/ropes_and_cables.md",
                "Swing" => "examples/swing.md",
                "Bodies in space" => "examples/space.md",
                "Gyroscopic effects" => "examples/gyroscopic_effects.md",
                "Wheels" => "examples/wheel.md",
                "Suspension systems" => "examples/suspension.md",
                "Quadrotor with cable-suspended load" => "examples/quad.md",
            ],
            "Components" => [
            "Frames" => "frames.md",
            "Joints" => "joints.md",
            "Components" => "components.md",
            "Forces" => "forces.md",
            "Sensors" => "sensors.md",
            "Trajectory planning" => "trajectory_planning.md",
            "Interfaces" => "interfaces.md",
            ],
            "Rotations and orientation" => "rotations.md",
            "3D rendering" => "rendering.md",
            "URDF import" => "urdf.md",
         ])

deploydocs(;
           repo = "github.com/JuliaComputing/Multibody.jl",
           devbranch = "main")
