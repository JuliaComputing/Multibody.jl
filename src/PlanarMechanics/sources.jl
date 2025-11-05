#=
within PlanarMechanics.Sources;
model WorldForce
  "External force and torque acting at frame_b, defined by 3 input signals and resolved in world frame"

  outer PlanarWorld planarWorld "Planar world model";

  parameter Modelica.Mechanics.MultiBody.Types.ResolveInFrameB resolveInFrame=
    Modelica.Mechanics.MultiBody.Types.ResolveInFrameB.frame_b
    "Frame in which output vector r_rel shall be resolved (1: world, 2: frame_b, 3: frame_resolve)";
  parameter Boolean animation=true "= true, if animation shall be enabled";

  parameter Real N_to_m(unit="N/m") = planarWorld.defaultN_to_m
    "Force arrow scaling (length = force/N_to_m)"
    annotation (Dialog(tab="Animation",group="If animation = true", enable=animation));
  parameter Real Nm_to_m(unit="N.m/m") = planarWorld.defaultNm_to_m
    "Torque arrow scaling (length = torque/Nm_to_m)"
    annotation (Dialog(tab="Animation",group="If animation = true", enable=animation));

  input SI.Diameter diameter = planarWorld.defaultArrowDiameter
    "Has no longer an effect and is only kept for backwards compatibility (arrow visualization by Vector now)"
    annotation (Dialog(tab="Animation",group="If animation = true", enable=animation));
  parameter SI.Length zPosition = planarWorld.defaultZPosition
    "Position z of cylinder representing the fixed translation"
    annotation (Dialog(tab="Animation",group="If animation = true", enable=animation));
  input MB.Types.Color color= PlanarMechanics.Types.Defaults.ForceColor
    "Color of arrow"
    annotation (HideResult=true, Dialog(tab="Animation",group="If animation = true",colorSelector=true,  enable=animation));
  input MB.Types.SpecularCoefficient specularCoefficient = planarWorld.defaultSpecularCoefficient
    "Reflection of ambient light (= 0: light is completely absorbed)"
    annotation (HideResult=true, Dialog(tab="Animation",group="If animation = true",enable=animation));

  Interfaces.Frame_b frame_b
    "Coordinate system fixed to the component with one cut-force and cut-torque"
    annotation (Placement(transformation(extent={{84,-16},{116,16}})));

  Modelica.Blocks.Interfaces.RealInput force[3]
    "x-, y-coordinates of force and torque resolved in world frame"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));

  Interfaces.Frame_resolve frame_resolve(fx = 0, fy = 0, t = 0, phi = phi) if resolveInFrame == Modelica.Mechanics.MultiBody.Types.ResolveInFrameB.frame_resolve
    "Coordinate system in which vector is optionally resolved, if useExtraFrame is true"
    annotation (Placement(transformation(extent={{-16,-16},{16,16}},rotation=90,origin={0,-100})));

  Real R[2,2] "Rotation matrix from world frame to frame_b";
  SI.Angle phi "Rotation angle of the additional frame_c";

protected
  SI.Position f_in_m[3]={force[1],force[2],0}/N_to_m "Force mapped from N to m for animation";
  SI.Position t_in_m[3]={0,0,force[3]}/Nm_to_m "Torque mapped from N.m to m for animation";
  MB.Frames.Orientation R_0 = MB.Frames.absoluteRotation(planarWorld.R,MB.Frames.planarRotation({0, 0, 1},phi,der(phi)));
  SI.Position rvisobj[3] = MB.Frames.resolve1(planarWorld.R,{frame_b.x,frame_b.y,zPosition}) + planarWorld.r_0;

  MB.Visualizers.Advanced.Vector arrowForce(
    coordinates=f_in_m,
    color=color,
    specularCoefficient=specularCoefficient,
    r=rvisobj,
    quantity=Modelica.Mechanics.MultiBody.Types.VectorQuantity.Force,
    headAtOrigin=true,
    R=R_0) if planarWorld.enableAnimation and animation;
  MB.Visualizers.Advanced.Vector arrowTorque(
    coordinates=t_in_m,
    color=color,
    specularCoefficient=specularCoefficient,
    r=rvisobj,
    quantity=Modelica.Mechanics.MultiBody.Types.VectorQuantity.Torque,
    headAtOrigin=true,
    R=planarWorld.R) if planarWorld.enableAnimation and animation;

equation
  if resolveInFrame == Modelica.Mechanics.MultiBody.Types.ResolveInFrameB.frame_b then
    phi = frame_b.phi;
  elseif resolveInFrame == Modelica.Mechanics.MultiBody.Types.ResolveInFrameB.world then
    phi = 0;
  end if;

  R = {{cos(phi), -sin(phi)}, {sin(phi),cos(phi)}};
  {frame_b.fx,frame_b.fy} + R*{force[1], force[2]} = {0, 0};
  frame_b.t + force[3]= 0;

  annotation (
    Icon(
      coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
      graphics={
        Line(
          points={{0,0},{0,-100}},
          color={95,95,95},
          pattern=LinePattern.Dot,
          visible=resolveInFrame==Modelica.Mechanics.MultiBody.Types.ResolveInFrameB.frame_resolve),
        Polygon(
          points={{-100,10},{20,10},{20,41},{90,0},{20,-41},{20,-10},{-100,-10},
              {-100,10}},
          lineColor={0,127,0},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-150,80},{150,40}},
          textString="%name",
          textColor={0,0,255})}),
    Documentation(
      info="<html>
<p>
The <strong>3</strong> signals of the <strong>force</strong> connector contain force and torque.
The first and second signal are interpreted as the x- and y-coordinates of
a&nbsp;<strong>force</strong> and the third is torque, acting at the frame connector
to which <code>frame_b</code> of this component is attached.
Note that torque is a&nbsp;scalar quantity, which is exerted perpendicular
to the x-y plane.
</p>
<p>
An example of this model is given in the following figure:
</p>

<div>
<img src=\"modelica://PlanarMechanics/Resources/Images/Sources/WorldForce.png\" alt=\"Modelica diagram\">
</div>

<p>
The parameter resolveInFrame defines in which frame the input force shall be resolved.
</p>
</html>",
      revisions="<html>
<p>
<img src=\"modelica://PlanarMechanics/Resources/Images/dlr_logo.png\" alt=\"DLR logo\">
<strong>Developed 2010 at the DLR Institute of System Dynamics and Control</strong>
</p>
</html>"));
end WorldForce;
=#

@component function WorldForceTorque(; name, resolve_frame = :world, render = true,
                                      scale = 0.1, color = [0, 1, 0, 0.5],
                                      radius = 0.05, z = 0)
    pars = @parameters begin
        render = render, [description = "Render the force in animations"]
        scale = scale, [description = "Scale factor for force visualization"]
        color[1:4] = color, [description = "Color of the force in animations"]
        radius = radius, [description = "Radius of the force arrow in animations"]
        z = z, [description = "z coordinate of the force in animations"]
    end

    systems = @named begin
        frame_b = Frame()
        force = RealInput(; nin = 2)
        torque = RealInput()
    end

    vars = @variables begin
        phi(t), [state_priority = -1, description = "Angle of resolution frame"]
    end

    # Calculations from begin block with conditional logic
    if resolve_frame === :world
        R = 1
    else
        R = ori_2d(frame_b)
    end
    # R = [cos(phi) -sin(phi); sin(phi) cos(phi)]

    equations = [
        if resolve_frame === :world
            phi ~ 0
        else
            phi ~ frame_b.phi
        end
        [frame_b.fx, frame_b.fy] + R * [force.u[1], force.u[2]] ~ [0, 0]
        frame_b.tau + torque.u ~ 0
    ]

    return System(equations, t; name, systems)
end