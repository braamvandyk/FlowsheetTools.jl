using MacroTools, Optim, JSON3, StructTypes, CSV, DataFrames

# Dicts for periodic table with atomic number => symbol => atomic mass
include("Atoms.jl")
using .Atoms

include("prettyround.jl")
include("components.jl")
include("streams.jl")
include("unitops.jl")
include("boundaries.jl")
include("kpis.jl")
include("closure.jl")

syscomps = Component[]
sysstreams = Stream[]
sysunitOps = UnitOp[]
sysboundaries = BalanceBoundary[]

