using MacroTools, Optim, JSON3, StructTypes, CSV, DataFrames, Statistics, Dates
import Base.setindex!
import Base.getindex
import Base.length



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


