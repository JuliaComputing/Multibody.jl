function World(; name)
    @named frame_b = Frame(name = name)
    @parameters n[1:3]=[0, 1, 0] g=9.81
    eqs = Equation[frame_b.r_0 .~ 0
                   vec(frame_b.R .~ nullrotation())]
    compose(ODESystem(eqs, t, [], [n; g]; name), frame_b)
end
