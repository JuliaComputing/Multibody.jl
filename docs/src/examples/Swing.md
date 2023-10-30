# Swing



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
        spring = Spring(c = inv(0.04/60))
        damper = Damper(d = 10.0)
        # spring = SpringDamperParallel(c = inv(0.04/60), d = 1, s_unstretched=0)
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

@mtkmodel Swing begin # Simple
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
@named model = Swing()
# ssys = structural_simplify(model, allow_parameters=false)
# prob = ODEProblem(ssys, ModelingToolkit.missing_variable_defaults(ssys), (0, 1))
# sol = solve(prob, Rodas4())



ssys = structural_simplify(IRSystem(model))
prob = ODEProblem(ssys, [], (0, 200))
sol = solve(prob, Rodas4(), abstol=1e-8, reltol=1e-8)
plot(sol, layout=21, size=(1900, 1000), legend=false, link=:x)


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
        body  = Body(m=6, isroot=true, I_11=1, I_22=1, I_33=1)

        damper = Damper(d=50.0)
    end
    @equations begin
        connect(world.frame_b, upper_trans1.frame_a, upper_trans2.frame_a)
        connect(rope1.frame_a, rope2.frame_a, upper_trans1.frame_b)
        connect(rope3.frame_a, rope4.frame_a, upper_trans2.frame_b)
        connect(rope1.frame_b, rope2.frame_b, rope3.frame_b, rope4.frame_b, body.frame_a)

        connect(world.frame_b, damper.frame_a)
        connect(body.frame_a, damper.frame_b)
    end
end

@named model = Swing()
model = complete(model)
ssys = structural_simplify(IRSystem(model))

prob = ODEProblem(ssys, [
    collect(model.body.r_0) .=> [0, -2, -0.5];
], (0, 10))
prob.u0[2] = -2
prob.u0[1] = -0.5


sol = solve(prob, Rodas4())
@assert SciMLBase.successful_retcode(sol)

plot(sol, idxs = [body.r_0...])
```


## 3D animation
Multibody.jl supports automatic 3D rendering of mechanisms, we use this feature to illustrate the result of the simulation below:

```@example spring_mass_system
import CairoMakie
Multibody.render(model, sol; z = -5, filename = "swing.gif") # Use "swing.mp4" for a video file
nothing # hide
```

![animation](swing.gif)