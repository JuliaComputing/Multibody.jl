using ModelingToolkit
using Multibody
using OrdinaryDiffEq # Contains the ODE solver we will use
using Plots
t = Multibody.t


JULIASIM_PURPLE = [87,87,219,255.0f0]./255 # RGBA
length_scale = 0.5 # This controls the frequency of the oscillations, smaller scale -> faster oscillations
radius_small = length_scale*0.2
radius_large = length_scale*0.3

@component function Logo(; name)
    pars = @parameters begin
    end

    systems = @named begin
        world = World()
        revl  = Revolute(; radius = radius_large, color=JULIASIM_PURPLE, axisflange=true)
        revl2 = Revolute(; radius = radius_large, color=JULIASIM_PURPLE, axisflange=true)
        revr  = Revolute(; radius = radius_small, color=JULIASIM_PURPLE, axisflange=true)
        bodyl = Body(m=1, radius = radius_small, color=JULIASIM_PURPLE)
        bodyr = Body(m=1, radius = radius_large, color=JULIASIM_PURPLE)
        bar_top = FixedTranslation(r=length_scale*[1, 0.05, 0], radius=length_scale*0.025, color=JULIASIM_PURPLE)
        barl   = FixedTranslation(r=length_scale*[1, -1, 0], radius=length_scale*0.025, color=JULIASIM_PURPLE)
        barr   = FixedTranslation(r=1.1*length_scale*[-1, -1, 0], radius=length_scale*0.025, color=JULIASIM_PURPLE)

        damperl  = Rotational.Damper(d=0.1)
        damperl2 = Rotational.Damper(d=0.01)
        damperr  = Rotational.Damper(d=0.01)
    end

    vars = @variables begin
    end

    equations = Equation[
        connect(revl.frame_a, world.frame_b)

        connect(revl.frame_b, barl.frame_a)
        connect(barl.frame_b, bodyl.frame_a)

        connect(world.frame_b, revl2.frame_a)
        connect(revl2.frame_b, bar_top.frame_a)
        connect(bar_top.frame_b, revr.frame_a)
        connect(revr.frame_b, barr.frame_a)
        connect(barr.frame_b, bodyr.frame_a)

        connect(revl.axis, damperl.flange_a)
        connect(revl2.axis, damperl2.flange_a)
        connect(revr.axis, damperr.flange_a)

        connect(revl.support, damperl.flange_b)
        connect(revl2.support, damperl2.flange_b)
        connect(revr.support, damperr.flange_b)
    ]

    return System(equations, t; name, systems)
end

@named logo = Logo()
logo = complete(logo)
ssys = multibody(logo)
prob = ODEProblem(ssys, [], (0.0, 3.51))
sol = solve(prob, Rodas5P())
Plots.plot(sol)

import GLMakie
framerate = 30
timevec = [zeros(30); range(sol.t[1], sol.t[end], step=1/framerate)] |> reverse
render(logo, sol, timevec; z=-2, x=-0.3, y=0.3, filename="JuliaSim_logo.gif", framerate)