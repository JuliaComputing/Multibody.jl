"""
Library to model planar mechanical multi-body systems inspired by https://github.com/dzimmer/PlanarMechanics
"""

module PlanarMechanics

using LinearAlgebra
import ModelingToolkitStandardLibrary.Mechanical.Rotational
import ModelingToolkitStandardLibrary.Mechanical.TranslationalModelica
import ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkit
using ...Blocks: RealInput, RealOutput
import ...@symcheck
import ..Multibody

export Frame, FrameResolve, PartialTwoFrames, ZeroPosition, ori_2d
include("utils.jl")

export Fixed, Body, BodyShape, FixedTranslation, Spring, Damper, SpringDamper
export SlipBasedWheelJoint, SimpleWheel, IdealPlanetary, DifferentialGear, OneDOFSlippingWheelJoint, OneDOFRollingWheelJoint
include("components.jl")

export Revolute, Prismatic
include("joints.jl")

export AbsolutePosition,
       RelativePosition, AbsoluteVelocity, RelativeVelocity, AbsoluteAcceleration,
       RelativeAcceleration, connect_sensor
include("sensors.jl")
end
