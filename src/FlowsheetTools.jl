module FlowsheetTools

export  Component, ComponentList, @comp, writecomponent, readcomponent, readcomponentlist!,
        Stream, StreamHistory, StreamList, StreamHistoryList, @stream, copystream!, deletestream!, renamestream!,
        copystreamhistory!, deletestreamhistory!, renamestreamhistory!, readstreamhistory,
        UnitOp, UnitOpHistory, UnitOpList, UnitOpHistoryList,
        BalanceBoundary, BalanceBoundaryHistory,
        calccorrections, closemb,
        conversion, selectivity


using CSV, DataFrames, JSON3, MacroTools, Optim, StructTypes, Statistics, Dates
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

end
