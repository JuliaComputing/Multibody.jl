using ModelingToolkit
using InteractiveUtils
using ModelingToolkit: t_nounits as t, D_nounits as D

versioninfo()
@variables x(t)=0.0 u(t)[1:3]
@named unused = System(Equation[], t, [u], [])
@named sys = System([D(x) ~ 1.0], t; systems = [unused])
compiled = mtkcompile(sys)
ODEProblem(compiled, [], (0.0, 1.0))
