# within StewartPlatform.Components;

# model Disc "Cylindrical disc used to model the base and the platform"
#     extends StewartPlatform.Icons.Disc; // Icon
    
# // Imports
#     import StewartPlatform.Types.Units.Direction;

# // Parameters
#     parameter StewartPlatform.Types.DiscParameters discParameters(D=1, alpha=from_deg(30)) "Parameters of the disc";
#     parameter Direction direction = Direction.down "Cylinder direction from x-y plane along z-axis";

#     protected
#         parameter SI.Position J[6,3] = StewartPlatform.Functions.jointsCoordinates(discParameters) "Coordinates of the joints resolved in disc frame";

# // Models
#     public
#         Modelica.Mechanics.MultiBody.Parts.FixedRotation fixedRotation1(r={J[1, 1],J[1, 2],J[1, 3]}, animation=false) annotation (Placement(transformation(extent={{20,10},{40,30}})));
#         Modelica.Mechanics.MultiBody.Parts.FixedRotation fixedRotation2(r={J[2, 1],J[2, 2],J[2, 3]}, animation=false) annotation (Placement(transformation(extent={{-10,-10},{10,10}},rotation=90,origin={0,50})));
#         Modelica.Mechanics.MultiBody.Parts.FixedRotation fixedRotation3(r={J[3, 1],J[3, 2],J[3, 3]}, animation=false) annotation (Placement(transformation(extent={{-10,-10},{10,10}},rotation=90,origin={-40,50})));
#         Modelica.Mechanics.MultiBody.Parts.FixedRotation fixedRotation4(r={J[4, 1],J[4, 2],J[4, 3]}, animation=false) annotation (Placement(transformation(extent={{10,-10},{-10,10}},rotation=90,origin={-40,-50})));
#         Modelica.Mechanics.MultiBody.Parts.FixedRotation fixedRotation5(r={J[5, 1],J[5, 2],J[5, 3]}, animation=false) annotation (Placement(transformation(extent={{10,-10},{-10,10}},rotation=90,origin={0,-50})));
#         Modelica.Mechanics.MultiBody.Parts.FixedRotation fixedRotation6(r={J[6, 1],J[6, 2],J[6, 3]}, animation=false) annotation (Placement(transformation(extent={{20,-30},{40,-10}})));

#         Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a "Frame in the center of the disc" annotation (Placement(transformation(extent={{-16,-16},{16,16}},rotation=90,origin={0,-100}), iconTransformation(extent={{-16,-16},{16,16}},rotation=90,origin={0,-74})));
#         Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b[6] "Frames where the joints are located" annotation (Placement(transformation(extent={{-16,-16},{16,16}},rotation=90,origin={0,100}), iconTransformation(extent={{-16,-16},{16,16}},rotation=90,origin={0,74})));
        
#         Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(
#             r=if direction == Direction.up then {0,0,discParameters.thickness} else {0,0,-discParameters.thickness},
#             lengthDirection=if direction == Direction.up then {0,0,1} else {0,0,-1},
#             length=discParameters.thickness,
#             diameter=discParameters.De,
#             density=discParameters.density)
#         annotation (Placement(transformation(extent={{-20,-10},{-40,10}})));

#     equation
#         assert(discParameters.D<=discParameters.De,
#         "Disc: the disc external diameter 'De' should be bigger or equal to the 'D' diameter of the circumference where the joints are located.",
#         level = AssertionLevel.warning);
        
#     connect(fixedRotation1.frame_a, bodyCylinder.frame_a) annotation(Line(points = {{20, 20}, {20, 20}, {20, 0}, {-20, 0}, {-20, 0}}));
#     connect(fixedRotation2.frame_a, bodyCylinder.frame_a) annotation(Line(points = {{0, 40}, {0, 40}, {0, 0}, {-20, 0}, {-20, 0}}));
#     connect(fixedRotation3.frame_a, bodyCylinder.frame_a) annotation(Line(points = {{-40, 40}, {-20, 40}, {-20, 0}, {-20, 0}}));
#     connect(fixedRotation4.frame_a, bodyCylinder.frame_a) annotation(Line(points = {{-40, -40}, {-20, -40}, {-20, 0}, {-20, 0}}));
#     connect(fixedRotation5.frame_a, bodyCylinder.frame_a) annotation(Line(points = {{0, -40}, {0, -40}, {0, 0}, {-20, 0}, {-20, 0}}));
#     connect(fixedRotation6.frame_a, bodyCylinder.frame_a) annotation(Line(points = {{20, -20}, {20, -20}, {20, 0}, {-20, 0}, {-20, 0}}));
    
#     connect(fixedRotation1.frame_b, frame_b[1]) annotation(Line(points = {{40, 20}, {40, 20}, {40, 100}, {0, 100}, {0, 100}}));
#     connect(fixedRotation2.frame_b, frame_b[2]) annotation(Line(points = {{0, 60}, {0, 60}, {0, 100}, {0, 100}}));
#     connect(fixedRotation3.frame_b, frame_b[3]) annotation(Line(points = {{-40, 60}, {-40, 60}, {-40, 100}, {0, 100}, {0, 100}}));
#     connect(fixedRotation4.frame_b, frame_b[4]) annotation(Line(points = {{-40, -60}, {-60, -60}, {-60, 100}, {0, 100}, {0, 100}}));
#     connect(fixedRotation5.frame_b, frame_b[5]) annotation(Line(points = {{0, -60}, {80, -60}, {80, 100}, {0, 100}, {0, 100}}));
#     connect(fixedRotation6.frame_b, frame_b[6]) annotation(Line(points = {{40, -20}, {60, -20}, {60, 100}, {0, 100}, {0, 100}}));
    
#     connect(frame_a, bodyCylinder.frame_a) annotation(Line(points = {{0, -100}, {-20, -100}, {-20, 0}, {-20, 0}}));
# end Disc;


@mtkmodel Disc begin # "Cylindrical disc used to model the base and the platform"
    
    # Parameters
    @parameters begin
        discParameters = DiscParameters(D=1, alpha=deg2rad(30)) "Parameters of the disc"
        J[1:6,1:3] = StewartPlatform.Functions.jointsCoordinates(discParameters) "Coordinates of the joints resolved in disc frame"
    end

    # Models
    @components begin
        fixedRotation1 = FixedRotation(r=[J[1, 1],J[1, 2],J[1, 3]], render=false)
        fixedRotation2 = FixedRotation(r=[J[2, 1],J[2, 2],J[2, 3]], render=false)
        fixedRotation3 = FixedRotation(r=[J[3, 1],J[3, 2],J[3, 3]], render=false)
        fixedRotation4 = FixedRotation(r=[J[4, 1],J[4, 2],J[4, 3]], render=false)
        fixedRotation5 = FixedRotation(r=[J[5, 1],J[5, 2],J[5, 3]], render=false)
        fixedRotation6 = FixedRotation(r=[J[6, 1],J[6, 2],J[6, 3]], render=false)

        frame_a #"Frame in the center of the disc"
        [frame_b[i] for i = 1:6] #"Frames where the joints are located"
        
        bodyCylinder = BodyCylinder(
            r = [0,0,discParameters.thickness]
            lengthDirection = [0,0,1],
            length = discParameters.thickness,
            diameter = discParameters.De,
            density = discParameters.density)
    end

    @equations begin
        # assert(discParameters.D<=discParameters.De,
        # "Disc: the disc external diameter 'De' should be bigger or equal to the 'D' diameter of the circumference where the joints are located.",
        # level = AssertionLevel.warning);
        
        connect(fixedRotation1.frame_a, bodyCylinder.frame_a)
        connect(fixedRotation2.frame_a, bodyCylinder.frame_a)
        connect(fixedRotation3.frame_a, bodyCylinder.frame_a)
        connect(fixedRotation4.frame_a, bodyCylinder.frame_a)
        connect(fixedRotation5.frame_a, bodyCylinder.frame_a)
        connect(fixedRotation6.frame_a, bodyCylinder.frame_a)
        
        connect(fixedRotation1.frame_b, frame_b[1])
        connect(fixedRotation2.frame_b, frame_b[2])
        connect(fixedRotation3.frame_b, frame_b[3])
        connect(fixedRotation4.frame_b, frame_b[4])
        connect(fixedRotation5.frame_b, frame_b[5])
        connect(fixedRotation6.frame_b, frame_b[6])
        
        connect(frame_a, bodyCylinder.frame_a)
    end
end 

# within StewartPlatform.Components;

# model Base "Base of a Stewart Platform"
#     extends StewartPlatform.Icons.Base; // Icon
#     extends Disc(final direction=StewartPlatform.Types.Units.Direction.down, final discParameters=if useGlobalParameters then gp.base else base);

# // Parameters
#     outer StewartPlatform.Components.GlobalParameters gp;

#     parameter Boolean useGlobalParameters = true "=true, if you want to use the base parameters defined in the 'gp' (GlobalParameters) object" annotation (choices(checkBox=true));
#     parameter StewartPlatform.Types.DiscParameters base(D=0.74, alpha=from_deg(12)) "Parameters of the base";
#     parameter SI.Position basePos[3] = {0,0,0} "Vector from world frame to base.frame_a resolved in world frame";
  
# // Models
#     Modelica.Mechanics.MultiBody.Parts.Fixed fixed(animation=false) annotation (Placement(transformation(extent={{-10,-11},{10,11}},rotation=90,origin={-80,-9})));        
#     Modelica.Mechanics.MultiBody.Parts.FixedRotation fixedRotation(animation=false, r=if useGlobalParameters then gp.basePos else basePos) annotation (Placement(transformation(extent={{-52,8},{-32,28}})));

# equation
#   connect(fixed.frame_b, fixedRotation.frame_a) annotation(Line(points = {{-80, 2}, {-80, 2}, {-80, 18}, {-52, 18}, {-52, 18}}));
#   connect(fixedRotation.frame_b, bodyCylinder.frame_a) annotation(Line(points = {{-32, 18}, {-20, 18}, {-20, 0}, {-20, 0}}, color = {95, 95, 95}));

# end Base;

@mtkmodel Base begin
    @parameters begin
        basePos[1:3] = [0,0,0] "Vector from world frame to base.frame_a resolved in world frame"
        gp = GlobalScope(GlobalParameters)
    end

    @extends bodyCylinder = disc = Disc(discParameters=gp.base)

    @components begin
        fixed = Fixed(render=false)
        fixedRotation = FixedRotation(r=basePos, render=false)
    end

    @equations begin
        connect(fixed.frame_b, fixedRotation.frame_a)
        connect(fixedRotation.frame_b, bodyCylinder.frame_a)
    end
end

# model Platform "Platform of a Stewart Platform"
#   extends StewartPlatform.Icons.Platform; // Icon
#   extends Disc(final direction=StewartPlatform.Types.Units.Direction.up, final discParameters = if useGlobalParameters then gp.platform else platform,
#       bodyCylinder(
#       final r_0(each fixed=fix_initPlatformPos, start = if useGlobalParameters then gp.initPlatformPos+gp.basePos else initPlatformPos+basePos),
#       final v_0(each fixed=fix_initPlatformVel),
#       final a_0(each fixed=fix_initPlatformAcc),
#       final angles_fixed=fix_initPlatformOrientation,
#       final angles_start=initAngles,
#       final sequence_start=initSequence,
#       final w_0_fixed=fix_initPlatformAngularVel,
#       final w_0_start=initPlatformAngularVel,
#       final z_0_fixed=fix_initPlatformAngularAcc,      
#       final z_0_start=initPlatformAngularAcc,
#       final enforceStates=_enforceStates,
#       final useQuaternions=_useQuaternions,
#       final sequence_angleStates=_sequence_angleStates));

# // Imports
#     import Modelica.Mechanics.MultiBody.Frames;
#     import TY = Modelica.Mechanics.MultiBody.Types;

# // Parameters
#     outer StewartPlatform.Components.GlobalParameters gp;

#     parameter Boolean useGlobalParameters = true "=true, if you want to use the parameters defined in the 'gp' (GlobalParameters) object"
#     parameter StewartPlatform.Types.DiscParameters platform(D=0.44, alpha=from_deg(100)) "Platform parameters";         
               
# // Initialization
#     parameter Boolean fix_initPlatformPos = false "=true, if you want to use the platform initial position as initial equation"
#     parameter Boolean fix_initPlatformVel = false "=true, if you want to use the platform initial velocity as initial equation"
#     parameter Boolean fix_initPlatformAcc = false "=true, if you want to use the platform initial acceleration as initial equation"
#     parameter Boolean fix_initPlatformOrientation = false "=true, if you want to use the platform initial orientation as initial equation"
#     parameter Boolean fix_initPlatformAngularVel = false "=true, if you want to use the platform initial angular velocity as initial equation"
#     parameter Boolean fix_initPlatformAngularAcc = false "=true, if you want to use the platform initial angular acceleration as initial equation"
#     parameter Boolean _enforceStates = false "= true, if absolute variables of body object shall be used as states (StateSelect.always)"
#     parameter Boolean _useQuaternions = false "= true, if quaternions shall be used as potential states otherwise use 3 angles as potential states"

#     parameter SI.Position basePos[3] = {0,0,0} "Vector from world frame to base.frame_a resolved in world frame"
#     parameter SI.Position initPlatformPos[3] = {0,0,1.45} "Coordinates of the platform frame resolved in base frame"
#     final parameter SI.Velocity initPlatformVel[3] = {0,0,0} "Initial platform velocity";
#     final parameter SI.Acceleration initPlatformAcc[3] = {0,0,0} "Initial platform acceleration";
#     final parameter SI.AngularVelocity initPlatformAngularVel[3]={0,0,0} "Initial platform angular velocity";
#     final parameter SI.AngularAcceleration initPlatformAngularAcc[3]={0,0,0} "Initial platform angular acceleration";
#     final parameter Modelica.Mechanics.MultiBody.Types.RotationSequence _sequence_angleStates={1,2,3} "Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states";
    
#     parameter TY.RotationTypes rotationType=TY.RotationTypes.RotationAxis "Type of rotation description"
    
#     parameter TY.Axis n={1,0,0} "Axis of rotation in base frame (= same as in base frame and platform frame)"
#     parameter NonSI.Angle_deg angle=0 "Angle to rotate base frame around axis n into platform frame"

#     parameter TY.Axis n_x={1,0,0} "Vector along x-axis of platform frame resolved in base frame"
#     parameter TY.Axis n_y={0,1,0} "Vector along y-axis of platform frame resolved in base frame"

#     parameter TY.RotationSequence sequence(min={1,1,1},max={3,3,3}) = {1,2,3} "Sequence of rotations"
#     parameter NonSI.Angle_deg angles[3]={0,0,0} "Rotation angles around the axes defined in 'sequence'"

#     //Conversion of the orientation in a sequence of rotations
#     final parameter Frames.Orientation R_rel=if useGlobalParameters
#                                              then Frames.nullRotation()
#                                              else if rotationType == TY.RotationTypes.RotationAxis
#                                              then
#                                              Frames.planarRotation(Modelica.Math.Vectors.normalizeWithAssert(n),from_deg(angle),0)
#                                              else if rotationType == TY.RotationTypes.TwoAxesVectors
#                                              then
#                                              Frames.from_nxy(n_x, n_y) else Frames.axesRotations(sequence, from_deg(angles),zeros(3))
#                                              "Fixed rotation object from base frame to platform frame";

#     final parameter TY.RotationSequence initSequence=if useGlobalParameters then gp.initSequence else {1,2,3} "Sequence of rotations to rotate base frame into platform frame; used to initialize the model";

#     final parameter SI.Angle initAngles[3]=if useGlobalParameters then gp.initAngles else Frames.axesRotationsAngles(R_rel,initSequence) "Initial values of angles to rotate base frame around 'initSequence' axes into platform frame; used to initialize the model";

# end Platform;

@mtkmodel Platform begin
    @parameters begin
        gp = GlobalScope(GlobalParameters)
        useGlobalParameters = true, [description = "If true, the parameters defined in the 'gp' (GlobalParameters) object are used"]
        platform = DiscParameters(D=0.44, alpha=deg2rad(100)), [description = "Platform parameters"]
        fix_initPlatformPos = false, [description = "If true, the platform initial position is used as initial equation"]
        fix_initPlatformVel = false, [description = "If true, the platform initial velocity is used as initial equation"]
        fix_initPlatformAcc = false, [description = "If true, the platform initial acceleration is used as initial equation"]
        fix_initPlatformOrientation = false, [description = "If true, the platform initial orientation is used as initial equation"]
        fix_initPlatformAngularVel = false, [description = "If true, the platform initial angular velocity is used as initial equation"]
        fix_initPlatformAngularAcc = false, [description = "If true, the platform initial angular acceleration is used as initial equation"]
        _enforceStates = false, [description = "If true, absolute variables of body object are used as states (StateSelect.always)"]
        _useQuaternions = false, [description = "If true, quaternions are used as potential states otherwise use 3 angles as potential states"]
        basePos[1:3] = [0,0,0], [description = "Vector from world frame to base.frame_a resolved in world frame"]
        initPlatformPos[1:3] = [0,0,1.45], [description = "Coordinates of the platform frame resolved in base frame"]
        initPlatformVel[1:3] = [0,0,0], [description = "Initial platform velocity"]
        initPlatformAcc[1:3] = [0,0,0], [description = "Initial platform acceleration"]
        initPlatformAngularVel[1:3] = [0,0,0], [description = "Initial platform angular velocity"]
        initPlatformAngularAcc[1:3] = [0,0,0], [description = "Initial platform angular acceleration"]
        _sequence_angleStates = [1,2,3], [description = "Sequence of rotations to rotate world frame into frame_a around the 3 angles used as potential states"]
        rotationType = RotationTypes.RotationAxis, [description = "Type of rotation description"]
        n = [1,0,0], [description = "Axis of rotation in base frame (= same as in base frame and platform frame)"]
        angle = 0, [description = "Angle to rotate base frame around axis n into platform frame"]
        n_x = [1,0,0], [description = "Vector along x-axis of platform frame resolved in base frame"]
        n_y = [0,1,0], [description = "Vector along y-axis of platform frame resolved in base frame"]
        sequence = [1,1,1], [description = "Sequence of rotations"]
        angles = [0,0,0], [description = "Rotation angles around the axes defined in 'sequence'"]
        R_rel = null_rotation(), [description = "Fixed rotation object from base frame to platform frame"]    
        initSequence = [1,2,3], [description = "Sequence of rotations to rotate base frame into platform frame; used to initialize the model"]
        initAngles = axes_rotations_angles(R_rel,initSequence), [description = "Initial values of angles to rotate base frame around 'initSequence' axes into platform frame; used to initialize the model"]
    end

    @extends bodyCylinder = disc = Disc(discParameters=platform)

    @components begin
        fixed = Fixed(render=false)
        fixedRotation = FixedRotation(r=basePos, render=false)
    end
end

# partial model PartialElectricCylinder "Partial model of an electric cylinder"

# // Parameters
#   parameter StewartPlatform.Types.ElectricCylinderParameters electricCylinderParameters "Parameters of the electric cylinder";  

# // Variables
#   SI.Length stroke "Current stroke, (=length-minPistonLength-boxLength-workingStroke/2)";
#   SI.Length length "Total length, distance between frame_a and frame_b";
#   SI.Velocity velocity "First derivative of 'length' (relative velocity)";
#   SI.Acceleration acceleration "Second derivative of 'length' (relative acceleration)";
#   SI.Force feedForce "Actuation force, in piston direction";
#   Real revolutions "Total number of rotations of the spindle for the given position (revolutions=0 <-> length=0)";
#   NonSI.AngularVelocity_rpm rotationalSpeed "Rotational speed of the spindle";

# // Models
#   Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (Placement(transformation(extent={{-116,-16},{-84,16}})));
#   Modelica.Mechanics.MultiBody.Joints.Prismatic 
#     prismatic(
#       stateSelect=electricCylinderParameters.stateSelect,
#       animation=false,
#       n={0,0,1},
#       useAxisFlange=true,
#       s(fixed = electricCylinderParameters.initialLengthFixed, start=electricCylinderParameters.initialLength),        
#       v(fixed = electricCylinderParameters.initialVelFixed, start=0),
#       a(fixed = electricCylinderParameters.initialAccFixed, start=0))
#     annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
#   Modelica.Mechanics.MultiBody.Parts.BodyBox 
#     bodyBox(
#       r={0,0,electricCylinderParameters.boxLength},
#       width=electricCylinderParameters.boxWidth,
#       density=(electricCylinderParameters.totalMass - electricCylinderParameters.movingMass)/(electricCylinderParameters.boxWidth^2*electricCylinderParameters.boxLength)) 
#     annotation (Placement(transformation(extent={{-60,-34},{-40,-14}})));
#   Modelica.Mechanics.MultiBody.Parts.BodyCylinder 
#     bodyCylinder(
#       r={0,0,electricCylinderParameters.pistonTotLength},
#       diameter=electricCylinderParameters.pistonDiameter,
#       density=4*electricCylinderParameters.movingMass/(electricCylinderParameters.pistonDiameter^2*pi*electricCylinderParameters.pistonTotLength),
#       color={255,0,0})
#     annotation (Placement(transformation(extent={{40,-34},{60,-14}})));
#   Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b annotation (Placement(transformation(extent={{84,-16},{116,16}})));

# equation
#   assert(abs(stroke)<=electricCylinderParameters.workingStroke/2, "PartialElectricCylinder: The current stroke is over the working stroke.", level = AssertionLevel.warning);
#   assert(electricCylinderParameters.minLength<electricCylinderParameters.maxLength,"PartialElectricCylinder: Limits must be consistent; however minLength>=maxLength");
#   assert(electricCylinderParameters.initialLength<=electricCylinderParameters.maxLength and electricCylinderParameters.initialLength>=electricCylinderParameters.minLength, "PartialElectricCylinder: Initial length is not within [minLength, maxLength]",level = AssertionLevel.warning);

#   length = prismatic.s;
#   velocity = prismatic.v;
#   acceleration = prismatic.a;
#   feedForce = prismatic.f;
#   stroke=length-electricCylinderParameters.boxLength-electricCylinderParameters.minPistonLength-electricCylinderParameters.workingStroke/2;
#   revolutions=length*1000/electricCylinderParameters.spindlePitch;
#   rotationalSpeed=velocity*1000/electricCylinderParameters.spindlePitch;

#   connect(frame_a, prismatic.frame_a) annotation (Line(points={{-100,0},{-10,0}},color={95,95,95},thickness=0.5));
#   connect(prismatic.frame_b, frame_b) annotation (Line(points={{10,0},{100,0}},color={95,95,95},thickness=0.5));
#   connect(bodyBox.frame_a, prismatic.frame_a) annotation (Line(points={{-60,-24},{-60,0},{-10,0}},color={95,95,95},thickness=0.5));
#   connect(bodyCylinder.frame_b, frame_b) annotation (Line(points={{60,-24},{60,0},{100,0}},color={95,95,95},thickness=0.5));
  
# end PartialElectricCylinder;

@mtkmodel PartialElectricCylinder begin
    @parameters begin
        electricCylinderParameters = ElectricCylinderParameters() "Parameters of the electric cylinder"
    end

    @variables begin
        stroke(t), [description = "Current stroke, (=length-minPistonLength-boxLength-workingStroke/2)"]
        length(t), [description = "Total length, distance between frame_a and frame_b"]
        velocity(t), [description = "First derivative of 'length' (relative velocity)"]
        acceleration(t), [description = "Second derivative of 'length' (relative acceleration)"]
        feedForce(t), [description = "Actuation force, in piston direction"]
        revolutions(t), [description = "Total number of rotations of the spindle for the given position (revolutions=0 <-> length=0)"]
        rotationalSpeed(t), [description = "Rotational speed of the spindle"]
    end

    @extends begin
        prismatic = Prismatic(
            # stateSelect=electricCylinderParameters.stateSelect,
            render=false,
            n=[0,0,1],
            axisflange=true,
            s=electricCylinderParameters.initialLength,
            v=0,
            a=0)
        bodyBox = BodyBox(
            r=[0,0,electricCylinderParameters.boxLength],
            width=electricCylinderParameters.boxWidth,
            density=(electricCylinderParameters.totalMass - electricCylinderParameters.movingMass)/(electricCylinderParameters.boxWidth^2*electricCylinderParameters.boxLength))
        bodyCylinder = BodyCylinder(
            r=[0,0,electricCylinderParameters.pistonTotLength],
            diameter=electricCylinderParameters.pistonDiameter,
            density=4*electricCylinderParameters.movingMass/(electricCylinderParameters.pistonDiameter^2*pi*electricCylinderParameters.pistonTotLength),
            color=[255,0,0,1])
    end

    @equations begin
        # assert(abs(stroke)<=electricCylinderParameters.workingStroke/2, "PartialElectricCylinder: The current stroke is over the working stroke.", level = AssertionLevel.warning)
        # assert(electricCylinderParameters.minLength<electricCylinderParameters.maxLength,"PartialElectricCylinder: Limits must be consistent; however minLength>=maxLength")
        # assert(electricCylinderParameters.initialLength<=electricCylinderParameters.maxLength and electricCylinderParameters.initialLength>=electricCylinderParameters.minLength, "PartialElectricCylinder: Initial length is not within [minLength, maxLength]",level = AssertionLevel.warning)

        length ~ prismatic.s
        velocity ~ prismatic.v
        acceleration ~ prismatic.a
        feedForce ~ prismatic.f
        stroke ~ length - electricCylinderParameters.boxLength - electricCylinderParameters.minPistonLength - electricCylinderParameters.workingStroke/2
        revolutions ~ length*1000/electricCylinderParameters.spindlePitch
        rotationalSpeed ~ velocity*1000/electricCylinderParameters.spindlePitch

        connect(prismatic.frame_a, frame_a)
        connect(prismatic.frame_b, frame_b)
        connect(bodyBox.frame_a, prismatic.frame_a)
        connect(bodyCylinder.frame_b, frame_b)
    end
end

# model ElectricCylinder "Mechanical linear drive unit with piston rod"
#   extends StewartPlatform.Icons.ElectricCylinder; // Icon
#   extends StewartPlatform.Components.PartialElectricCylinder;

# // Models  
#   SpindleDrive spindleDrive(spindlePitch=electricCylinderParameters.spindlePitch) annotation (Placement(transformation(extent={{-18,26},{2,46}})));
#   Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a annotation (Placement(transformation(extent={{-110,56},{-90,76}}), iconTransformation(extent={{-110,56},{-90,76}})));
#   Modelica.Mechanics.Rotational.Components.IdealGear idealGear(ratio=electricCylinderParameters.ratio) annotation (Placement(transformation(extent={{-48,26},{-28,46}})));

# equation

#   connect(spindleDrive.flangeT, prismatic.axis) annotation (Line(points={{2,36},{8,36},{8,6}}, color={0,127,0}));
#   connect(idealGear.flange_b, spindleDrive.flangeR) annotation (Line(points={{-28,36},{-18,36}}, color={0,0,0}));
#   connect(flange_a, idealGear.flange_a) annotation (Line(points={{-100,66},{-100,36},{-48,36}}, color={0,0,0}));

# end ElectricCylinder;

@mtkmodel ElectricCylinder begin
    @extends PartialElectricCylinder

    @components begin
        spindleDrive = SpindleDrive(spindlePitch=electricCylinderParameters.spindlePitch)
        flange_a = Flange()
        idealGear = IdealGear(ratio=electricCylinderParameters.ratio)
    end

    @equations begin
        connect(spindleDrive.flangeT, prismatic.axis)
        connect(idealGear.flange_b, spindleDrive.flangeR)
        connect(flange_a, idealGear.flange_a)
    end

end

# model SpindleDrive "Spindle drive transforming rotational into translational motion"
#   extends StewartPlatform.Icons.SpindleDrive; // Icon

# // Imports
#   import StewartPlatform.Types.Units.Pitch;

# // Parameters
#   parameter Pitch spindlePitch = 5 "Spindle pitch";
#   final parameter Real ratio = 1000 * 2 * pi / spindlePitch "Transmission ratio (flangeR.phi/flangeT.s)";

# // Models
#   Modelica.Mechanics.Rotational.Interfaces.Flange_a flangeR "Flange of rotational shaft" annotation(Placement(transformation(extent = {{-110, -10}, {-90, 10}})));
#   Modelica.Mechanics.Translational.Interfaces.Flange_b flangeT "Flange of translational rod" annotation(Placement(transformation(extent = {{90, -10}, {110, 10}})));
#   Modelica.Mechanics.Rotational.Components.IdealGearR2T rotationToTranslation(ratio = ratio) annotation(Placement(visible = true, transformation(extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
# equation
#   connect(flangeR, rotationToTranslation.flangeR) annotation(Line(points = {{-100, 0}, {-10, 0}, {-10, 0}, {-10, 0}}));
#   connect(rotationToTranslation.flangeT, flangeT) annotation(Line(points = {{10, 0}, {98, 0}, {98, 0}, {100, 0}}, color = {0, 127, 0}));

# end SpindleDrive;

@mtkmodel SpindleDrive begin
    @parameters begin
        spindlePitch = 5, [description = "Spindle pitch"]
        ratio = 1000 * 2 * pi / spindlePitch, [description = "Transmission ratio (flangeR.phi/flangeT.s)"]
    end

    @components begin
        flangeR = Flange() 
        flangeT = Flange()
        rotationToTranslation = IdealGearR2T(ratio=ratio)
    end

    @equations begin
        connect(flangeR, rotationToTranslation.flangeR)
        connect(rotationToTranslation.flangeT, flangeT)
    end
end


# model Leg "Six-degree-of-freedom leg with electric cylinder and servo motor"
#   extends StewartPlatform.Icons.Leg; // Icon

# // Imports
#   import StewartPlatform.Types.*;

# // Parameters
#     // Universal joint
#     parameter UniversalJointParameters universalJointParameters "Parameters of the universal joint"

#     // Electric Cylinder
#     parameter ElectricCylinderParameters electricCylinderParameters "Parameters of the electric cylinder"

#     // Servo Motor
#     parameter ServoMotorParameters servoMotorParameters "Parameters of the servo motor"

#     // Spherical joint
#     parameter SphericalJointParameters sphericalJointParameters "Parameters of the spherical joint"

# // Models
#   Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a "Frame for connection with a base's joint frame"
#   Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b "Frame for connection with a platform's joint frame"

#   ElectricCylinder electricCylinder(electricCylinderParameters=electricCylinderParameters)

#   UniversalJoint universalJoint(universalJointParameters = universalJointParameters)

#   ServoMotor servoMotor(servoMotorParameters = servoMotorParameters)

#   SphericalJoint sphericalJoint(sphericalJointParameters=sphericalJointParameters)
  
#   StewartPlatform.Interfaces.AxisControlBus axisControlBus
    
# equation
#   connect(frame_a, universalJoint.frame_a)
#   connect(universalJoint.frame_b, electricCylinder.frame_a)
#   connect(electricCylinder.frame_b, sphericalJoint.frame_a)
#   connect(sphericalJoint.frame_b, frame_b)
#   connect(servoMotor.flange, electricCylinder.flange_a)
#   connect(axisControlBus, servoMotor.axisControlBus)
# end Leg;

@mtkmodel Leg begin
    @parameters begin
        universalJointParameters = UniversalJointParameters(), [description="Parameters of the universal joint"]
        electricCylinderParameters = ElectricCylinderParameters(), [description="Parameters of the electric cylinder"]
        servoMotorParameters = ServoMotorParameters(), [description="Parameters of the servo motor"]
        sphericalJointParameters = SphericalJointParameters(), [description="Parameters of the spherical joint"]
    end

    @components begin
        frame_a = Frame()
        frame_b = Frame()
        electricCylinder = ElectricCylinder(electricCylinderParameters=electricCylinderParameters)
        universalJoint = UniversalJoint(universalJointParameters = universalJointParameters)
        servoMotor = ServoMotor(servoMotorParameters = servoMotorParameters)
        sphericalJoint = SphericalJoint(sphericalJointParameters=sphericalJointParameters)
        axisControlBus = AxisControlBus()
    end

    @equations begin
        connect(frame_a, universalJoint.frame_a)
        connect(universalJoint.frame_b, electricCylinder.frame_a)
        connect(electricCylinder.frame_b, sphericalJoint.frame_a)
        connect(sphericalJoint.frame_b, frame_b)
        connect(servoMotor.flange, electricCylinder.flange_a)
        connect(axisControlBus, servoMotor.axisControlBus)
    end
end

# model ServoMotor "Simplified model of a permanent magnet synchronous induction machine"
#   extends StewartPlatform.Icons.ServoMotor; // Icon

#   // Imports
#   import Modelica.Constants.inf;

#   // Parameters
#   parameter StewartPlatform.Types.ServoMotorParameters servoMotorParameters "Parameters of the servo motor";

#   // Variables
#   SI.Angle phi "Absolute mechanical rotation angle (flange.phi)";
#   SI.AngularVelocity w "Mechanical angular velocity (der(phi))";
#   SI.AngularAcceleration a "Mechanical angular acceleration (der(w))";
#   SI.Torque Tref "Reference mechanical torque for shaft";
#   SI.Torque T "Mechanical torque for shaft (flange.tau)";
#   SI.Power Pmecc "Output mechanical power";

#   // Models
#   Modelica.Mechanics.Rotational.Interfaces.Flange_b flange "Shaft"
#   StewartPlatform.Interfaces.AxisControlBus axisControlBus

# protected
#   Modelica.Blocks.Nonlinear.Limiter limiter(uMax = if servoMotorParameters.enableLimiter then servoMotorParameters.Tmax else inf, uMin = if servoMotorParameters.enableLimiter then -servoMotorParameters.Tmax else -inf) 
#   Modelica.Mechanics.Rotational.Sources.Torque torque
#   Modelica.Mechanics.Rotational.Sensors.AngleSensor angleSensor
#   Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor
#   Modelica.Mechanics.Rotational.Sensors.AccSensor accSensor

# equation
#   phi = angleSensor.phi;
#   w = speedSensor.w;
#   a = accSensor.a;
#   Tref = limiter.u;
#   T = limiter.y;
#   Pmecc = w * T;
  
#   connect(torque.flange, flange)
#   connect(angleSensor.flange, flange)
#   connect(speedSensor.flange, flange)
#   connect(accSensor.flange, flange)
#   connect(limiter.y, torque.tau)
#   connect(axisControlBus.refTorque, limiter.u)
#   connect(angleSensor.phi, axisControlBus.angularPos)
#   connect(speedSensor.w, axisControlBus.angularVel)
#   connect(accSensor.a, axisControlBus.angularAcc)
#   connect(limiter.y, axisControlBus.appliedTorque);
  
# end ServoMotor;


@mtkmodel ServoMotor begin
    @parameters begin
        servoMotorParameters = ServoMotorParameters(), [description="Parameters of the servo motor"]
    end

    @variables begin
        phi(t), [description="Absolute mechanical rotation angle (flange.phi)"]
        w(t), [description="Mechanical angular velocity (der(phi))"]
        a(t), [description="Mechanical angular acceleration (der(w))"]
        Tref(t), [description="Reference mechanical torque for shaft"]
        T(t), [description="Mechanical torque for shaft (flange.tau)"]
        Pmecc(t), [description="Output mechanical power"]
    end

    @components begin
        flange = Flange()
        axisControlBus = AxisControlBus()
        limiter = Limiter(uMax = if servoMotorParameters.enableLimiter then servoMotorParameters.Tmax else Inf, uMin = if servoMotorParameters.enableLimiter then -servoMotorParameters.Tmax else -Inf)
        torque = Torque()
        angleSensor = AngleSensor()
        speedSensor = SpeedSensor()
        accSensor = AccSensor()
    end

    @equations begin
        phi ~ angleSensor.phi
        w ~ speedSensor.w
        a ~ accSensor.a
        Tref ~ limiter.u
        T ~ limiter.y
        Pmecc ~ w * T

        connect(torque.flange, flange)
        connect(angleSensor.flange, flange)
        connect(speedSensor.flange, flange)
        connect(accSensor.flange, flange)
        connect(limiter.y, torque.tau)
        connect(axisControlBus.refTorque, limiter.u)
        connect(angleSensor.phi, axisControlBus.angularPos)
        connect(speedSensor.w, axisControlBus.angularVel)
        connect(accSensor.a, axisControlBus.angularAcc)
        connect(limiter.y, axisControlBus.appliedTorque)
    end
end


# model SixLegs "Set of six legs of a Stewart Platform appropiately connected, composed of electric cylinders and servo motors"
#     extends StewartPlatform.Icons.SixLegs; // Icon
# // Parameters
#     outer StewartPlatform.Components.GlobalParameters gp;

#     // Universal joints
#     parameter UniversalJointParameters universalJointParameters[6] = gp.universalJointParameters "Parameters of the universal joints"

#     // Electric Cylinders
#     parameter ElectricCylinderParameters electricCylinderParameters[6] = gp.electricCylinderParameters "Parameters of the electric cylinders"

#     // Servo Motors
#     parameter ServoMotorParameters servoMotorParameters[6] = gp.servoMotorParameters "Parameters of the servo motors"

#     // Spherical joints
#     parameter SphericalJointParameters sphericalJointParameters[6] = gp.sphericalJointParameters "Parameters of the spherical joints"

# // Variables
#     SI.Torque motorTorqueDIS[6] "Display torque applied by motors";

# // Models
#   Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_base[6] "Frames for connection with the base's joints"
#   Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_platform[6] "Frames for connection with the platform's joints"
#   Interfaces.ControlBus controlBus
  
#   Leg leg1
#     (
#       universalJointParameters = universalJointParameters[1],
#       electricCylinderParameters = electricCylinderParameters[1],
#       servoMotorParameters = servoMotorParameters[1],
#       sphericalJointParameters = sphericalJointParameters[1]
#     )
       
#   Leg leg2
#     (
#       universalJointParameters = universalJointParameters[2],
#       electricCylinderParameters = electricCylinderParameters[2],
#       servoMotorParameters = servoMotorParameters[2],
#       sphericalJointParameters = sphericalJointParameters[2]
#     )

#   Leg leg3
#     (
#       universalJointParameters = universalJointParameters[3],
#       electricCylinderParameters = electricCylinderParameters[3],
#       servoMotorParameters = servoMotorParameters[3],
#       sphericalJointParameters = sphericalJointParameters[3]
#     )
   
#   Leg leg4
#     (
#       universalJointParameters = universalJointParameters[4],
#       electricCylinderParameters = electricCylinderParameters[4],
#       servoMotorParameters = servoMotorParameters[4],
#       sphericalJointParameters = sphericalJointParameters[4]
#     )

#   Leg leg5
#     (
#       universalJointParameters = universalJointParameters[5],
#       electricCylinderParameters = electricCylinderParameters[5],
#       servoMotorParameters = servoMotorParameters[5],
#       sphericalJointParameters = sphericalJointParameters[5]
#     )
   
#   Leg leg6
#     (
#       universalJointParameters = universalJointParameters[6],
#       electricCylinderParameters = electricCylinderParameters[6],
#       servoMotorParameters = servoMotorParameters[6],
#       sphericalJointParameters = sphericalJointParameters[6]
#     )
   

# equation
#   motorTorqueDIS[1]=leg1.servoMotor.T;
#   motorTorqueDIS[2]=leg2.servoMotor.T;
#   motorTorqueDIS[3]=leg3.servoMotor.T;
#   motorTorqueDIS[4]=leg4.servoMotor.T;
#   motorTorqueDIS[5]=leg5.servoMotor.T;
#   motorTorqueDIS[6]=leg6.servoMotor.T;

#   connect(frame_platform[1], leg1.frame_b)
#   connect(frame_platform[2], leg2.frame_b)
#   connect(frame_platform[3], leg3.frame_b)
#   connect(frame_platform[4], leg4.frame_b)
#   connect(frame_platform[5], leg5.frame_b)
#   connect(frame_platform[6], leg6.frame_b)

#   connect(frame_base[1], leg1.frame_a)
#   connect(frame_base[2], leg2.frame_a)
#   connect(frame_base[3], leg3.frame_a)
#   connect(frame_base[4], leg4.frame_a)
#   connect(frame_base[5], leg5.frame_a)
#   connect(frame_base[6], leg6.frame_a)

#   connect(controlBus.axisControlBus1, leg1.axisControlBus)
#   connect(controlBus.axisControlBus2, leg2.axisControlBus)
#   connect(controlBus.axisControlBus3, leg3.axisControlBus)
#   connect(controlBus.axisControlBus4, leg4.axisControlBus)
#   connect(controlBus.axisControlBus5, leg5.axisControlBus)
#   connect(controlBus.axisControlBus6, leg6.axisControlBus)
    
# end SixLegs;

@mtkmodel SixLegs begin
    @parameters begin
        gp = GlobalScope(GlobalParameters)
        universalJointParameters[1:6] = gp.universalJointParameters, [description="Parameters of the universal joints"]
        electricCylinderParameters[1:6] = gp.electricCylinderParameters, [description="Parameters of the electric cylinders"]
        servoMotorParameters[1:6] = gp.servoMotorParameters, [description="Parameters of the servo motors"]
        sphericalJointParameters[1:6] = gp.sphericalJointParameters, [description="Parameters of the spherical joints"]
    end

    @variables begin
        motorTorqueDIS(t)[1:6], [description="Display torque applied by motors"]
    end

    @components begin
        frame_base[1:6] = Frame()
        frame_platform[1:6] = Frame()
        controlBus = ControlBus()
        leg1 = Leg(
            universalJointParameters = universalJointParameters[1],
            electricCylinderParameters = electricCylinderParameters[1],
            servoMotorParameters = servoMotorParameters[1],
            sphericalJointParameters = sphericalJointParameters[1]
        )
        leg2 = Leg(
            universalJointParameters = universalJointParameters[2],
            electricCylinderParameters = electricCylinderParameters[2],
            servoMotorParameters = servoMotorParameters[2],
            sphericalJointParameters = sphericalJointParameters[2]
        )
        leg3 = Leg(
            universalJointParameters = universalJointParameters[3],
            electricCylinderParameters = electricCylinderParameters[3],
            servoMotorParameters = servoMotorParameters[3],
            sphericalJointParameters = sphericalJointParameters[3]
        )
        leg4 = Leg(
            universalJointParameters = universalJointParameters[4],
            electricCylinderParameters = electricCylinderParameters[4],
            servoMotorParameters = servoMotorParameters[4],
            sphericalJointParameters = sphericalJointParameters[4]
        )
        leg5 = Leg(
            universalJointParameters = universalJointParameters[5],
            electricCylinderParameters = electricCylinderParameters[5],
            servoMotorParameters = servoMotorParameters[5],
            sphericalJointParameters = sphericalJointParameters[5]
        )
        leg6 = Leg(
            universalJointParameters = universalJointParameters[6],
            electricCylinderParameters = electricCylinderParameters[6],
            servoMotorParameters = servoMotorParameters[6],
            sphericalJointParameters = sphericalJointParameters[6]
        )
    end
    @equations begin
        motorTorqueDIS[1] ~ leg1.servoMotor.T
        motorTorqueDIS[2] ~ leg2.servoMotor.T
        motorTorqueDIS[3] ~ leg3.servoMotor.T
        motorTorqueDIS[4] ~ leg4.servoMotor.T
        motorTorqueDIS[5] ~ leg5.servoMotor.T
        motorTorqueDIS[6] ~ leg6.servoMotor.T

        connect(frame_platform[1], leg1.frame_b)
        connect(frame_platform[2], leg2.frame_b)
        connect(frame_platform[3], leg3.frame_b)
        connect(frame_platform[4], leg4.frame_b)
        connect(frame_platform[5], leg5.frame_b)
        connect(frame_platform[6], leg6.frame_b)

        connect(frame_base[1], leg1.frame_a)
        connect(frame_base[2], leg2.frame_a)
        connect(frame_base[3], leg3.frame_a)
        connect(frame_base[4], leg4.frame_a)
        connect(frame_base[5], leg5.frame_a)
        connect(frame_base[6], leg6.frame_a)

        connect(controlBus.axisControlBus1, leg1.axisControlBus)
        connect(controlBus.axisControlBus2, leg2.axisControlBus)
        connect(controlBus.axisControlBus3, leg3.axisControlBus)
        connect(controlBus.axisControlBus4, leg4.axisControlBus)
        connect(controlBus.axisControlBus5, leg5.axisControlBus)
        connect(controlBus.axisControlBus6, leg6.axisControlBus)
    end
end

# model Controller "Cascade of controllers to control the axes"
# // Parameters
#     outer StewartPlatform.Components.GlobalParameters gp;
#     parameter Boolean useGlobalParameters = true "=true, if you want to use the base parameters defined in the 'gp' (GlobalParameters) object"
#     parameter StewartPlatform.Types.DiscParameters base(D=0.74, alpha=from_deg(12)) "Parameters of the base";
#     parameter StewartPlatform.Types.DiscParameters platform(D=0.44, alpha=from_deg(100)) "Parameters of the platform";
#     parameter StewartPlatform.Types.Units.Pitch spindlePitch = gp.electricCylinderParameters[1].spindlePitch "Spindle pitch of the electric cylinder for rotation-dispacement conversion";
#     parameter Real ratio = gp.electricCylinderParameters[1].ratio "Transmission ratio (servomotor.phi/spindle.phi)";
#     parameter SI.Torque initOutputs = 0.01107 "Initial value of the torque commands (Controller outputs)";

# // Controller
#     parameter Modelica.Blocks.Types.SimpleController controllerType = Modelica.Blocks.Types.SimpleController.PID "Type of controllers"
#     parameter Real P(min=0, unit="1") = gp.servoMotorParameters[1].Tmax/500 "Proportional action"
#     parameter SI.Time Ti(min=Modelica.Constants.small) = 0.5 "Time constant of Integrator blocks"
#     parameter SI.Time Td(min=0)= 0.01 "Time constant of Derivative blocks"
#     parameter SI.Time preFilterTimeConstant = 0.1 "Time constant of the prefilter applied to legth ref."
#     parameter SI.Time postFilterTimeConstant = 0.001 "Time constant of the postfilter (additional pole)"
#     parameter SI.Torque maxOutput = gp.servoMotorParameters[1].Tmax "The controllers output are limited within [-maxOutput,maxOutput]"

# // Inverse kinematic
#     parameter Boolean limitOutputs = true "=true, if you want to limit the outputs within [minLength,maxLength]"
#     parameter Boolean stopWhenSaturated = true "When a saturation is detected all output signals maintein their current values until all output signals return within the limits"
#     parameter SI.Length maxLength = gp.maxLength "Max leg length"
#     parameter SI.Length minLength = gp.minLength "Min leg length"
#     parameter SI.Time T_velocity(min=Modelica.Constants.small) = 0.01 "Time constants for velocity derivative (T>0 required; T=0 is ideal derivative block)"
#     parameter SI.Time T_acc(min=Modelica.Constants.small) = 0.1 "Time constants for acceleration derivative (T>0 required; T=0 is ideal derivative block)"
     

# // Models
#   Interfaces.Pose inputPose "Desired pose for the platform resolved in base frame"
#   ReferenceSignals.InverseKinematic inverseKinematic(
#     platform = if useGlobalParameters then gp.platform else platform,
#     base = if useGlobalParameters then gp.base else base,
#     limitOutputs=limitOutputs,
#     stopWhenSaturated=stopWhenSaturated,
#     maxLength=maxLength,
#     minLength=minLength,
#     T_velocity=T_velocity,
#     T_acc=T_acc)
       

#   Modelica.Blocks.Continuous.LimPID PID[6](
#     each yMax=maxOutput,
#     each controllerType=controllerType,
#     each k=P,
#     each Ti=Ti,
#     each Td=Td,
#     each y_start=initOutputs,
#     each initType=if controllerType <> Modelica.Blocks.Types.SimpleController.P then Modelica.Blocks.Types.Init.InitialOutput else Modelica.Blocks.Types.Init.NoInit)

#   Interfaces.ControlBus controlBus()

#   Modelica.Blocks.Math.Gain gain[6](each k=ratio*2*pi*1000/spindlePitch)

#   Modelica.Blocks.Continuous.FirstOrder preFilter[6](each T=preFilterTimeConstant, each initType=Modelica.Blocks.Types.Init.SteadyState)
#   Modelica.Blocks.Continuous.FirstOrder postFilter[6](each initType=Modelica.Blocks.Types.Init.SteadyState, each T=postFilterTimeConstant)
# equation
#   connect(inverseKinematic.pose, inputPose)
#   connect(inverseKinematic.legLength, preFilter.u)
#   connect(preFilter.y, gain.u)
#   connect(gain.y, PID.u_s)
#   connect(PID.y, postFilter.u)
  
#   connect(postFilter[1].y, controlBus.axisControlBus1.refTorque)
#   connect(postFilter[2].y, controlBus.axisControlBus2.refTorque);
#   connect(postFilter[3].y, controlBus.axisControlBus3.refTorque);
#   connect(postFilter[4].y, controlBus.axisControlBus4.refTorque);
#   connect(postFilter[5].y, controlBus.axisControlBus5.refTorque);
#   connect(postFilter[6].y, controlBus.axisControlBus6.refTorque);

#   connect(PID[1].u_m, controlBus.axisControlBus1.angularPos)
#   connect(PID[2].u_m, controlBus.axisControlBus2.angularPos);
#   connect(PID[3].u_m, controlBus.axisControlBus3.angularPos);
#   connect(PID[4].u_m, controlBus.axisControlBus4.angularPos);
#   connect(PID[5].u_m, controlBus.axisControlBus5.angularPos);
#   connect(PID[6].u_m, controlBus.axisControlBus6.angularPos);
  
# end Controller;

@mtkmodel Controller begin
    @parameters begin
        gp = GlobalScope(GlobalParameters)
        useGlobalParameters = true, [description="=true, if you want to use the base parameters defined in the 'gp' (GlobalParameters) object"]
        base = Disc(D=0.74, alpha=from_deg(12)), [description="Parameters of the base"]
        platform = Disc(D=0.44, alpha=from_deg(100)), [description="Parameters of the platform"]
        spindlePitch = gp.electricCylinderParameters[1].spindlePitch, [description="Spindle pitch of the electric cylinder for rotation-dispacement conversion"]
        ratio = gp.electricCylinderParameters[1].ratio, [description="Transmission ratio (servomotor.phi/spindle.phi)"]
        initOutputs = 0.01107, [description="Initial value of the torque commands (Controller outputs)"]
        controllerType = SimpleController.PID, [description="Type of controllers"]
        P = gp.servoMotorParameters[1].Tmax/500, [description="Proportional action"]
        Ti = 0.5, [description="Time constant of Integrator blocks"]
        Td = 0.01, [description="Time constant of Derivative blocks"]
        preFilterTimeConstant = 0.1, [description="Time constant of the prefilter applied to legth ref."]
        postFilterTimeConstant = 0.001, [description="Time constant of the postfilter (additional pole)"]
        maxOutput = gp.servoMotorParameters[1].Tmax, [description="The controllers output are limited within [-maxOutput,maxOutput]"]
        limitOutputs = true, [description="=true, if you want to limit the outputs within [minLength,maxLength]"]
        stopWhenSaturated = true, [description="When a saturation is detected all output signals maintein their current values until all output signals return within the limits"]
        maxLength = gp.maxLength, [description="Max leg length"]
        minLength = gp.minLength, [description="Min leg length"]
        T_velocity = 0.01, [description="Time constants for velocity derivative (T>0 required; T=0 is ideal derivative block)"]
        T_acc = 0.1, [description="Time constants for acceleration derivative (T>0 required; T=0 is ideal derivative block)"]
    end

    @components begin
        inputPose = Pose()
        inverseKinematic = InverseKinematic(
            platform = ifelse(useGlobalParameters, gp.platform, platform),
            base = ifelse(useGlobalParameters, gp.base, base),
            limitOutputs=limitOutputs,
            stopWhenSaturated=stopWhenSaturated,
            maxLength=maxLength,
            minLength=minLength,
            T_velocity=T_velocity,
            T_acc=T_acc)
        PID = LimPID(
            each yMax=maxOutput,
            each controllerType=controllerType,
            each k=P,
            each Ti=Ti,
            each Td=Td,
            each y_start=initOutputs,
            each initType=if controllerType <> SimpleController.P then Init.InitialOutput else Init.NoInit)
        controlBus = ControlBus(Axisbus = StewartAxisControlBus)
        gain = Gain(each k=ratio*2*pi*1000/spindlePitch)
        preFilter = FirstOrder(each T=preFilterTimeConstant, each initType=Init.SteadyState)
        postFilter = FirstOrder(each initType=Init.SteadyState, each T=postFilterTimeConstant)
    end

    @equations begin
        connect(inverseKinematic.pose, inputPose)
        connect(inverseKinematic.legLength, preFilter.u)
        connect(preFilter.y, gain.u)
        connect(gain.y, PID.u_s)
        connect(PID.y, postFilter.u)
        
        connect(postFilter[1].y, controlBus.axisControlBus1.refTorque)
        connect(postFilter[2].y, controlBus.axisControlBus2.refTorque)
        connect(postFilter[3].y, controlBus.axisControlBus3.refTorque)
        connect(postFilter[4].y, controlBus.axisControlBus4.refTorque)
        connect(postFilter[5].y, controlBus.axisControlBus5.refTorque)
        connect(postFilter[6].y, controlBus.axisControlBus6.refTorque)

        connect(PID[1].u_m, controlBus.axisControlBus1.angularPos)
        connect(PID[2].u_m, controlBus.axisControlBus2.angularPos)
        connect(PID[3].u_m, controlBus.axisControlBus3.angularPos)
        connect(PID[4].u_m, controlBus.axisControlBus4.angularPos)
        connect(PID[5].u_m, controlBus.axisControlBus5.angularPos)
        connect(PID[6].u_m, controlBus.axisControlBus6.angularPos)
    end

end


# model InverseKinematic "The outputs are the six leg lengths required to have the input pose"
#     extends StewartPlatform.Icons.InverseKinematic; # Icon

# # Imports
#     import ModelicaServices.Machine.eps; #this constant is used to verify if there is a saturation

# # Parameters
#     outer StewartPlatform.Components.GlobalParameters gp "Object with all global parameters";

#     parameter StewartPlatform.Types.DiscParameters platform = gp.platform "Parameters of the platform";
#     parameter StewartPlatform.Types.DiscParameters base = gp.base "Parameters of the base";

#     parameter Boolean limitOutputs = true "=true, if you want to limit the outputs within [minLength,maxLength]" 
#     parameter Boolean stopWhenSaturated = true "When a saturation is detected all output signals maintain their current values until all output signals return within the limits" 
#     parameter SI.Length maxLength = gp.maxLength "Max leg length" 
#     parameter SI.Length minLength = gp.minLength "Min leg length" 

#     parameter SI.Time T_velocity(min=Modelica.Constants.small) = 0.01 "Time constants for velocity derivative (T>0 required; T=0 is ideal derivative block)" 
#     parameter SI.Time T_acc(min=Modelica.Constants.small) = 0.1 "Time constants for acceleration derivative (T>0 required; T=0 is ideal derivative block)" 

# # Variables
#     Boolean saturationFlag "=true, when there is a saturation";
#     SI.Velocity legVelocity[6] "Leg velocities, =der(legLength)";
#     SI.Acceleration legAcceleration[6] "Leg accelerations, =der(legVelocity)";

# protected
#     SI.Length lengthRef[6] "Output before saturation, the result of inverse kinematic";
#     SI.Length lengthSat[6] "Output after saturation";
#     SI.Length outputSample[6] "When a saturation is detected this variables save the last valid outputs";

# # Models
# public
#   Interfaces.Pose pose "Input pose" 
#   Modelica.Blocks.Interfaces.RealOutput legLength[6] "Leg lengths for the input pose"

# protected
#   Modelica.Blocks.Nonlinear.Limiter limiter[6](each uMax=maxLength, each uMin=minLength) 

#   Modelica.Blocks.Continuous.Derivative posDerivative[6](each T=T_velocity, each initType=Modelica.Blocks.Types.Init.InitialOutput) 
#   Modelica.Blocks.Continuous.Derivative velDerivative[6](each T=T_acc, each initType=Modelica.Blocks.Types.Init.InitialOutput) 

#   Interfaces.PoseDeMux deMux 


# initial equation
#   outputSample = lengthSat; #Initialize the valid outputs with the input signals after saturation.
#                             #If 'stopWhenSaturated' is true and, at the beginning of the simulation, there are saturaions then there may be discontinuities in the outputs of this model

# equation
#   assert(maxLength>=minLength,"InverseKinematic: Limits must be consistent. However minLength>maxLength.");

#   lengthRef = StewartPlatform.Functions.legsLength(base,platform,deMux.positionOut,deMux.sequenceOut,deMux.orientationOut); #resolve inverse kinematic

#   limiter.u=lengthRef; #apply saturation
#   lengthSat=limiter.y;

#   if lengthRef[1]-eps>lengthSat[1] or lengthRef[2]-eps>lengthSat[2] or lengthRef[3]-eps>lengthSat[3] or
#      lengthRef[4]-eps>lengthSat[4] or lengthRef[5]-eps>lengthSat[5] or lengthRef[6]-eps>lengthSat[6] or
#      lengthRef[1]+eps<lengthSat[1] or lengthRef[2]+eps<lengthSat[2] or lengthRef[3]+eps<lengthSat[3] or
#      lengthRef[4]+eps<lengthSat[4] or lengthRef[5]+eps<lengthSat[5] or lengthRef[6]+eps<lengthSat[6] then
#                                                                                                         #check if there was a saturation
#        saturationFlag=true;
#    else
#      saturationFlag=false;
#   end if;

#   when saturationFlag then #immidiatly after a saturation save the last output
#            outputSample=lengthSat;
#   end when;

#   if limitOutputs then

#     if stopWhenSaturated and saturationFlag then #Select the output based on the parameters and possible saturations
#         legLength=outputSample;
#     else
#         legLength=lengthSat;
#     end if;

#   else
#     legLength=lengthRef;
#   end if;

#   posDerivative.u=legLength; #compute the output derivatives with appropriate blocks
#   legVelocity=posDerivative.y;

#   velDerivative.u=legVelocity;
#   legAcceleration=velDerivative.y;

#   connect(deMux.poseIn, pose) 
  
# end InverseKinematic;

# @mtkmodel InverseKinematic begin
#     @parameters begin
#         gp = GlobalScope(GlobalParameters)
#         platform = gp.platform, [description="Parameters of the platform"]
#         base = gp.base, [description="Parameters of the base"]
#         limitOutputs = true, [description="=true, if you want to limit the outputs within [minLength,maxLength]"]
#         stopWhenSaturated = true, [description="When a saturation is detected all output signals maintain their current values until all output signals return within the limits"]
#         maxLength = gp.maxLength, [description="Max leg length"]
#         minLength = gp.minLength, [description="Min leg length"]
#         T_velocity = 0.01, [description="Time constants for velocity derivative (T>0 required; T=0 is ideal derivative block)"]
#         T_acc = 0.1, [description="Time constants for acceleration derivative (T>0 required; T=0 is ideal derivative block)"]
#     end

#     @variables begin
#         saturationFlag(t), [description="=true, when there is a saturation"]
#         legVelocity(t)[1:6], [description="Leg velocities, =der(legLength)"]
#         legAcceleration(t)[1:6], [description="Leg accelerations, =der(legVelocity)"]
#         lengthRef(t) = Length[1:6], [description="Output before saturation, the result of inverse kinematic"]
#         lengthSat(t) = Length[1:6], [description="Output after saturation"]
#         outputSample(t) = Length[1:6], [description="When a saturation is detected this variables save the last valid outputs"]
#         pose(t) = Pose(), [description="Input pose"]
#         legLength(t) = RealOutput[1:6], [description="Leg lengths for the input pose"]
#     end

#     @components begin
#         limiter = Limiter(each uMax=maxLength, each uMin=minLength)
#         posDerivative = Derivative(each T=T_velocity, each initType=Init.InitialOutput)
#         velDerivative = Derivative(each T=T_acc, each initType=Init.InitialOutput)
#         deMux = PoseDeMux()
#     end

#     # @initial_equations begin
#     #     outputSample ~ lengthSat
#     # end

#     @equations begin
#         assert(maxLength>=minLength,"InverseKinematic: Limits must be consistent. However minLength>maxLength.")

#         lengthRef ~ StewartPlatform.Functions.legsLength(base,platform,deMux.positionOut,deMux.sequenceOut,deMux.orientationOut)

#         limiter.u ~ lengthRef
#         lengthSat ~ limiter.y

#         when saturationFlag then
#             outputSample ~ lengthSat
#         end

#         if limitOutputs then
#             if stopWhenSaturated and saturationFlag then
#                 legLength ~ outputSample
#             else
#                 legLength ~ lengthSat
#             end
#         else
#             legLength ~ lengthRef
#         end

#         posDerivative.u ~ legLength
#         legVelocity ~ posDerivative.y

#         velDerivative.u ~ legVelocity
#         legAcceleration ~ velDerivative.y

#         connect(deMux.poseIn, pose)
#     end
# end