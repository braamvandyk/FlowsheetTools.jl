#----------------------------------------------------------------------------
#
#----Mixer-------------------------------------------------------------------
#
#----------------------------------------------------------------------------


"""

    function mixer!(fs, outlets, inlets, params)

Calculation for mixer UnitOps. Combines all feed streams into a single outlet stream.
    streams::StreamList (should be automatically assigned to the UnitOp upon creation)
    outlets::Vector{String} or equivalent
    inlets::Vector{String} or equivalent
    params can by of any type relevant to the calculations

Will throw if more than one outlet stream is defined. Values in `params` are ignored.

"""
function mixer!(streams, outlets, inlets, params)
    @argcheck length(outlets) == 1 "mixers can have only one outlet stream"

    tempstream = sum(streams[inlets])

    streams[outlets[1]] = renamestream(tempstream, outlets[1])
    
    return nothing
end