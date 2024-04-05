using ModelingToolkit
using Multibody, JuliaSimCompiler
using OrdinaryDiffEq # Contains the ODE solver we will use
using Plots
t = Multibody.t
W(args...; kwargs...) = Multibody.world


JULIASIM_PURPLE = [87,87,219,255.0f0]./255 # RGBA
length_scale = 0.5
radius_small = length_scale*0.2
radius_large = length_scale*0.3

@mtkmodel Logo begin
    @components begin
        world = World()
        revl = Revolute(; radius = radius_large, color=JULIASIM_PURPLE)
        revr = Revolute(; radius = radius_small, color=JULIASIM_PURPLE)
        bodyl = Body(m=1, radius = radius_small, color=JULIASIM_PURPLE)
        bodyr = Body(m=1, radius = radius_large, color=JULIASIM_PURPLE)
        bar_top = FixedTranslation(r=length_scale*[1, 0.05, 0], radius=length_scale*0.025, color=JULIASIM_PURPLE)
        bar_l   = FixedTranslation(r=length_scale*[1, -1, 0], radius=length_scale*0.025, color=JULIASIM_PURPLE)
        bar_r   = FixedTranslation(r=1.1*length_scale*[-1, -1, 0], radius=length_scale*0.025, color=JULIASIM_PURPLE)
    end
    @equations begin
        connect(revl.frame_a, world.frame_b)

        connect(revl.frame_b, bar_l.frame_a)
        connect(bar_l.frame_b, bodyl.frame_a)

        connect(world.frame_b, bar_top.frame_a)
        connect(bar_top.frame_b, revr.frame_a)
        connect(revr.frame_b, bar_r.frame_a)
        connect(bar_r.frame_b, bodyr.frame_a)
    end
end

@named logo = Logo()
logo = complete(logo)
ssys = structural_simplify(IRSystem(logo))
prob = ODEProblem(ssys, [], (0.0, 3.51))
sol = solve(prob, Rodas5P())
Plots.plot(sol)

using GLMakie
render(logo, sol; z=-2, x=-0.3, y=0.3, filename="JuliaSim_logo.gif")