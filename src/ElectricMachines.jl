module ElectricMachines

using CairoMakie
using ForwardDiff
using LaTeXStrings
using OrdinaryDiffEq, DiffEqCallbacks
using SymPy

include("models/dc.jl")
include("plotutils.jl")
include("trajectory.jl")

end # module
