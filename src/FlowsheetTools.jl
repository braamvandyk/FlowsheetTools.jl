# TODO Should macros use ~ instead of ->??
# TODO replace all @error and @assert with exceptions (where relevant)
# TODO Add precompile workload

module FlowsheetTools

export  Component, ComponentList, @comp, writecomponent, readcomponent, readcomponentlist!, names,
        Stream, StreamList, @stream, copystream!, deletestream!, renamestream!, renamestream, emptystream, readstreamhistory, refreshcomplist,
        UnitOp, UnitOpList, @unitop, mixer!, flowsplitter!, componentplitter!,
        BalanceBoundary, @boundary, showdata, 
        calccorrections, closemb,
        conversion, molar_selectivity,
        Flowsheet, addunitop!, setorder!, generateBFD


using Base64,                   # Used to encode strings to pass Mermaid diagram definition to server
      ChangePointDetection,     # Used to detect changepoints (step changes)
      Dates,                    # Used with all DateTime types
      DelimitedFiles,           # Used for reading stream histories
      Downloads,                # Used to pass Mermaid diagram definition to server
      HypothesisTests,          # Used for Augmented Dickey-Fuller test to see if data has a slop between changepoints
      Interpolations,           # Used for extrapolation to endpoints in cleaning up data
      Loess,                    # Used for smoothing and filling data sets
      MacroTools,               # Used for @comp, @stream etc
      Missings,                 # Used for handling missing data in data cleanup
      Optim,                    # Used to minimize objective function for mass balance reconciliations (replace with JuMP?)
      OrderedCollections,       # Used for OrderedDict
      PrettyTables,             # Used for pretty printing
      Statistics,               # Basic statistical functions
      TimeSeries                # Used for TimeArray in Streams


import Base.setindex!
import Base.getindex
import Base.length
import Base.copy
import Base.iterate
import Base.show
import Base.:+
import Base.:*
# TODO Add these
# import Base.==
# import Base.≈



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
