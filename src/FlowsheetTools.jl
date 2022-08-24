module FlowsheetTools

export  Component, ComponentList, @comp, writecomponent, readcomponent, readcomponentlist!,
        Stream, StreamHistory, StreamList, StreamHistoryList, @stream, copystream!, deletestream!, renamestream!,
        copystreamhistory!, deletestreamhistory!, renamestreamhistory!, readstreamhistory, showdata,
        UnitOp, UnitOpHistory, UnitOpList, UnitOpHistoryList,
        BalanceBoundary, BalanceBoundaryHistory, @boundary, @boundaryhist,
        calccorrections, closemb,
        conversion, selectivity,
        @unitop, @unitophist, mixer!,
        Flowsheet, addunitop!, setorder!


using JSON3, MacroTools, Optim, StructTypes, Statistics, Dates, DelimitedFiles, PrettyTables
import Base.setindex!
import Base.getindex
import Base.length
import Base.copy



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
include("flowsheets.jl")

end
