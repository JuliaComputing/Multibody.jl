# Swing

In this example, we will model a swing consisting of a rigid seat suspended in 4 ropes, mounted symmetrically in a ceiling. Each rope is modeled as a stiff rod with a small point mass at the center of gravity, terminated by a parallel spring-damper to model slight flexibility in the ropes. The ceiling mounting points are modeled as spherical joints, i.e., they do not transmit any torque in any direction. The rim of the seat is modeled as 4 rigid bodies configured in a square, as well as one point mass representing the load, located slightly below the rim assembly.

![animation](swing.gif)


We start by defining a single rope component and attach it to a body in order to verify that it's working.

```@example SWING
using Multibody
using ModelingToolkit
using Plots
using JuliaSimCompiler
using OrdinaryDiffEq

t = Multibody.t
D = Differential(t)
world = Multibody.world

W(args...; kwargs...) = Multibody.world

@mtkmodel Rope begin
    @components begin
        frame_a = Frame()
        frame_b = Frame()
        joint1 = Spherical(isroot=true, enforceState=true)
        rope = BodyShape(r=[0.0,-1,0], m=0.05, isroot=false, radius=0.01)
        spring = Spring(c = inv(0.04/60), m=0.01)
        damper = Damper(d = 10.0)
    end
    @equations begin
        connect(frame_a, joint1.frame_a)
        connect(joint1.frame_b, rope.frame_a)

        connect(rope.frame_b, spring.frame_a, damper.frame_a)
        connect(spring.frame_b, damper.frame_b, frame_b)

        # connect(rope.frame_b, frame_b)
    end
end

# @named rope = Rope(rope.r = [0, -2, 0])

@mtkmodel SimpleSwing begin # Simple
    @structural_parameters begin
        h = 2
        w = 0.4
    end
    @components begin
        world = W()
        upper_trans1 = FixedTranslation(r=[-w/2, 0, 0])
        rope1 = Rope(rope.r=[-w/2, h, -w/2])
        body  = Body(m=6, isroot=true, I_11=0.1, I_22=0.1, I_33=0.1)
        damper = Damper(d=10.0)
    end
    @equations begin
        connect(world.frame_b, upper_trans1.frame_a)
        connect(rope1.frame_a, upper_trans1.frame_b)
        # connect(world.frame_b, rope1.frame_a)
        connect(rope1.frame_b, body.frame_a)
        
        connect(world.frame_b, damper.frame_a)
        connect(body.frame_a, damper.frame_b)
    end
end
@named model = SimpleSwing()
ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, [], (0, 4))
sol = solve(prob, Rodas4(), abstol=1e-8, reltol=1e-8)
plot(sol, layout=21, size=(1900, 1000), legend=false, link=:x)
```
This makes for a rather interesting-looking springy pendulum!

```@example SWING
import CairoMakie
Multibody.render(model, sol, 0:0.01:4; z = -5, filename = "simple_swing.gif") # Use "swing.mp4" for a video file
nothing # hide
```

Next, we create the full swing assembly

```@example SWING
@mtkmodel Swing begin
    @structural_parameters begin
        h = 2
        w = 0.4
    end
    @components begin
        world = W()
        upper_trans1 = FixedTranslation(r=[-w/2, 0, 0])
        upper_trans2 = FixedTranslation(r=[ w/2, 0, 0])
        rope1 = Rope(rope.r=[-w/2, -h, -w/2])
        rope2 = Rope(rope.r=[-w/2, -h,  w/2])
        rope3 = Rope(rope.r=[ w/2, -h, -w/2])
        rope4 = Rope(rope.r=[ w/2, -h,  w/2])

        body_back  = BodyShape(m=0.1, r = [w, 0, 0])
        body_front = BodyShape(m=0.1, r = [w, 0, 0])
        body_left  = BodyShape(m=0.1, r = [0, 0, w])
        body_right = BodyShape(m=0.1, r = [0, 0, -w])


        body  = Body(m=6, isroot=true, r_cm = [w/2, -w/2, w/2])

        damper = Damper(d=0.5)
    end
    @equations begin
        # Rope assembly
        connect(world.frame_b, upper_trans1.frame_a, upper_trans2.frame_a)
        connect(rope1.frame_a, rope2.frame_a, upper_trans1.frame_b)
        connect(rope3.frame_a, rope4.frame_a, upper_trans2.frame_b)

        # Body assembly
        connect(body_back.frame_a, body_left.frame_a, rope1.frame_b)
        connect(body_left.frame_b, body_front.frame_a, rope2.frame_b)
        connect(body_front.frame_b, body_right.frame_a, rope4.frame_b)
        connect(body_right.frame_b, rope3.frame_b) # Don't close the rigid kinematic loop
        connect(body_back.frame_a, body.frame_a)

        # World damping (damps swing motion)
        connect(world.frame_b, damper.frame_a)
        connect(body.frame_a, damper.frame_b)
    end
end

@named model = Swing()
model = complete(model)
ssys = structural_simplify(IRSystem(model))

d = 10
prob = ODEProblem(ssys, [
    collect(model.body_left.body.r_0) .=> [0, -2, -0.5];
    model.rope1.damper.d => d;
    model.rope2.damper.d => d;
    model.rope3.damper.d => d;
    model.rope4.damper.d => d;
], (0.0, 0.23))
# prob.u0[1] = -0.5
# prob.u0[2] = -2


@time sol = solve(prob, Rosenbrock23(autodiff=false))
@assert SciMLBase.successful_retcode(sol)

Plots.plot(sol, idxs = [model.body.r_0...])
```

```@example spring_mass_system
import CairoMakie
Multibody.render(model, sol; z = -5, filename = "swing.gif") # Use "swing.mp4" for a video file
nothing # hide
```

![animation](swing.gif)