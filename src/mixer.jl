#----------------------------------------------------------------------------
#
#----Mixer-------------------------------------------------------------------
#
#----------------------------------------------------------------------------


"""
    mixer!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)

Calculation for mixer UnitOps. Combines all feed streams into a single outlet stream.
Will error if more than one outlet stream is defined. Values in `params` are ignored.
"""
function mixer!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)
    @argcheck length(outlets) == 1 "mixers can have only one outlet stream"

    name = streamlist[outlets[1]].name     
    tempstream = sum(streamlist[inlets])

    tempstream = renamestream(tempstream, name)
    streamlist[name] = tempstream
    
    return nothing
end