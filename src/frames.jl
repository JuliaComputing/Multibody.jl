# function Oriantation(; name, R=collect(1.0I(3)), w=zeros(3))
#     # T: Transformation matrix from world frame to local frame
#     # w: Absolute angular velocity of local frame, resolved in local frame

#     @variables R(t)[1:3, 1:3]=R [description="Orientation rotation matrix ∈ SO(3)"]
#     @variables w(t)[1:3]=w [description="angular velocity"]
#     R,w = collect.((R,w))

#     ODESystem(Equation[], t, [vec(R); w], [], name = name)
# end

@connector function Frame(; name, derived_w = true)
    @variables r_0(t)[1:3] [description = "Position vector directed from the origin of the world frame to the connector frame origin, resolved in world frame"]
    @variables f(t)[1:3] [connect = Flow, description = "Cut force resolved in connector frame"]
    @variables tau(t)[1:3] [connect = Flow, description = "Cut torque resolved in connector frame"]
    r_0, f, tau = collect.((r_0, f, tau))
    # R: Orientation object to rotate the world frame into the connector frame
    
    R = NumRotationMatrix(; name, derived_w)

    ODESystem(Equation[], t, [r_0; f; tau], []; name, metadata=Dict(:orientation => R, :frame => true))
end



