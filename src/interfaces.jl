function PartialTwoFrames(; name)
    @named frame_a = Frame()
    @named frame_b = Frame()
    compose(ODESystem([], t; name), frame_a, frame_b)
end