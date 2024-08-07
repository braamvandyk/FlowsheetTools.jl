#----------------------------------------------------------------------------
#
#----Definitions-------------------------------------------------------------
#
#----------------------------------------------------------------------------


struct UnitOp
    name::String

    streams::StreamList
    # Connected streams with history
    inlets::Vector{String}
    outlets::Vector{String}

    f!::Function # An optional function that will write to outlets from inlets and params
    params # Default paramaters can be added upon creation of the object

    # Inner constructor to validate input data
    function UnitOp(name, streamlist, inlets, outlets, f!, params)
        numdata = streamlist[inlets[1]].numdata
        timestamps = timestamp(streamlist[inlets[1]].massflows)
        
        for inletname in inlets
            inlet = streamlist[inletname]
            inlet_ts = timestamp(inlet.massflows)

            inlet.numdata != numdata && throw(DimensionMismatch("all streams must must have identical history lengths."))
            !all(inlet_ts .== timestamps) && throw(ArgumentError("all streams must must have identical timestamps."))
        end
        for outletname in outlets
            outlet = streamlist[outletname]
            outlet_ts = timestamp(outlet.massflows)

            outlet.numdata != numdata && throw(DimensionMismatch("all streams must must have identical history lengths."))
            !all(outlet_ts .== timestamps) && throw(ArgumentError("all streams must must have identical timestamps."))
        end

        new(name, streamlist, inlets, outlets, f!, params)
    end
end


function (u::UnitOp)(newparams = nothing)
    if isnothing(newparams)
        u.f!(u.streams, u.outlets, u.inlets, u.params)
    else
        u.f!(u.streams, u.outlets, u.inlets, newparams)
    end
end


struct UnitOpList
    list::OrderedDict{String, UnitOp}
end


#----------------------------------------------------------------------------
#
#----Constructors------------------------------------------------------------
#
#----------------------------------------------------------------------------

"""
    UnitOp(name::String, streamlist::StreamList, inlets::Vector{String}, outlets::Vector{String})
    UnitOp(name::String, streamlist::StreamList, inlets::Vector{String}, outlets::Vector{String}, f!::Function)
    UnitOp(name::String, streamlist::StreamList, inlets::Vector{String}, outlets::Vector{String}, f!::Function, params)

Constructor for a UnitOp. It is recommended to rather use the @unitop macro.
"""
function UnitOp(name::String, streamlist::StreamList, inlets::Vector{String}, outlets::Vector{String})
    return UnitOp(name, streamlist, inlets, outlets, passive, nothing)
end


function UnitOp(name::String, streamlist::StreamList, inlets::Vector{String}, outlets::Vector{String}, f!::Function)
    return UnitOp(name, streamlist, inlets, outlets, f!, nothing)
end


"""
    UnitOpList()

Constructor for an empty UnitOpList. Unit operations are added when created with the @unitop macro.
"""
function UnitOpList()
    l = OrderedDict{String, UnitOp}()
    return UnitOpList(l)
end


#----------------------------------------------------------------------------
#
#----Active UnitOps----------------------------------------------------------
#
#----------------------------------------------------------------------------

"""
    passive(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)

Default calculation for UnitOps. Simply a placeholder - does no calculation.
"""
@inline function passive(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)
    return nothing
end

include("mixer.jl")
include("splitters.jl")
include("reactors.jl")


#----------------------------------------------------------------------------
#
#----Base overloads----------------------------------------------------------
#
#----------------------------------------------------------------------------


function Base.setindex!(A::UnitOpList, X::UnitOp, idx::String)
    if length(A.list) == 0
        A.list[idx] = X
    else
        # Verify that all the unit ops reference the same StreamList.
        # Get the first item in the `list` field. As `list` is a `Dict`, this returns a `Pair`.
        # We get the value entry using the `second` field of the `Pair`, which returns a `UnitOp`,
        # of which we get the `streamlist` field.
        currentlist = first(A.list).second.streams
        X.streams != currentlist && throw(ArgumentError("all unit operations in UnitOpList must reference the same StreamList"))

        A.list[idx] = X
    end
end


function Base.getindex(A::UnitOpList, idx::String)
    return A.list[idx]
end


function Base.getindex(A::UnitOpList, idxs::Vector{String})
    res = UnitOp[]
    for idx in idxs
        push!(res, A.list[idx])
    end
    return res
end


function Base.length(A::UnitOpList)
    return length(A.list)
end


# Pretty printing for UnitOp objects
function Base.show(io::IO, u::UnitOp)
    println(io, "Unit Operation:  $(u.name)")
    println(io, "Feed streams:    ", [s for s in u.inlets])
    println(io)
    println(io, "Product streams: ", [s for s in u.outlets])
    
    strm = first(u.streams).second
    
    # If streams hold only one data point, don't print this
    if strm.numdata > 1
        ts = timestamp(strm.massflows)
        println(io)
        println(io, "Data length:\t$(strm.numdata)")
        println(io, "Data starts:\t$(ts[begin])")
        println(io, "Data ends:\t$(ts[end])") 
    end
end


function  Base.show(io::IO, uol::UnitOpList)
    println(io, "UnitOperation list:")
    if length(uol.list) > 0
        for unitop in uol.list
            println(io, "  ", unitop.first)
        end
        
        uo = first(uol.list).second
        sl = uo.streams
        inlets = uo.inlets
        strm = sl[inlets[1]]
        
        # If streams hold only one data point, don't print this
        if strm.numdata > 1
            println(io)
            println(io, "Data length:\t$(strm.numdata)")
            ts = timestamp(strm.massflows)
            println(io, "Data starts:\t$(ts[begin])")
           println(io, "Data ends:\t$(ts[end])")
        end
    else
        println(io, "Empty list")
    end
end


function Base.iterate(A::UnitOpList)
    return iterate(A.list)
end


function Base.iterate(A::UnitOpList, state)
    return iterate(A.list, state)
end


#----------------------------------------------------------------------------
#
#----Macros------------------------------------------------------------------
#
#----------------------------------------------------------------------------


"""
Defines a `UnitOp with` the specified name and connected streams.

@unitop begin
    inlets --> ["H2", "C2"]
    outlets --> ["Mixed"]
    calc --> mixer!
end "Mixer" fs

This will create a UnitOp named "Mixer", saved into fs.unitops["Mixer"].
The specified streams refer to entries in fs.streams::StreamList.

Two optional parameters may be specified, i.e. calc and params:

calc is a function, e.g. 

    function mixer!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)

params is an iterable passed to the function.

Together, these specify a calculation to perform to calculate the outlet(s) from the inlet(s), when
calling `fs.unitops["Mixer"]()`

Alternative parameter values may also be specified in the call, e.g.
    fs.unitops["Mixer"]([1, 2, 3])

"""
macro unitop(ex::Expr, name::String, fs::Symbol)      
    local inlets = String[]
    local outlets = String[]
    local func = nothing
    local params = nothing
    
    for line in ex.args
        match_comp = @capture(line, inlets --> [ins__])
        if match_comp
            for strm in ins
                push!(inlets, strm)
            end
        end
        
        match_comp = @capture(line, outlets --> [outs__])
        if match_comp
            for strm in outs
                push!(outlets, strm)
            end
        end

        match_comp = @capture(line, calc --> unitopfunc_)
        if match_comp
            func = unitopfunc
        end

        match_comp = @capture(line, params --> unitopparams_)
        if match_comp
            params = unitopparams
        end

    end

    if isnothing(func)
        return quote
            $(esc(fs)).unitops[$name] = UnitOp($name, $(esc(fs)).streams, $inlets, $outlets)
            addtorun!($(esc(fs)), $name)
        end
    else    
        return quote
            $(esc(fs)).unitops[$name] = UnitOp($name, $(esc(fs)).streams, $inlets, $outlets, $(esc(func)), $(esc(params)))
            addtorun!($(esc(fs)), $name)
        end
    end
end


#----------------------------------------------------------------------------
# 
#----Utilities---------------------------------------------------------------
# 
#----------------------------------------------------------------------------


"""

    function deleteunitop!(fs, name)

Delete a stream from the Flowsheet's StreamList.
"""
function deleteunitop!(fs, name)
    @argcheck fs isa Flowsheet "fs must be a Flowsheet"
    if name in keys(fs.unitops.list)
        delete!(fs.unitops.list, name)

        idx = findfirst(==(name), fs.rununits)
        if !isnothing(idx)
            deleteat!(fs.rununits, idx)
            
            gap = fs.runorder[idx]
            deleteat!(fs.runorder, idx)
            
            # Fix gap in the run order
            for i in eachindex(fs.runorder)
                if fs.runorder[i] > gap
                    fs.runorder[i] -= 1
                end
            end
        end
        
        # Also clear the boundaries to ensure consistency
        
        for (bname, boundary) in fs.boundaries.list
            if name in boundary.included_units
                delete!(fs.boundaries.list, bname)
            end
        end
    end
    
    return nothing
end


"""

function deleteunitops!(fs)

Delete all streams from the Flowsheet's StreamList.

"""
function deleteunitops!(fs)
    @argcheck fs isa Flowsheet "fs must be a Flowsheet"
    for name in keys(fs.unitops.list)
        delete!(fs.unitops.list, name)
    end
    for i in eachindex(fs.rununits)
        deleteat!(fs.rununits, 1)
    end
    for i in eachindex(fs.runorder)
        deleteat!(fs.runorder, 1)
    end

    # Also clear the boundaries to ensure consistency

    for bname in keys(fs.boundaries.list)
        delete!(fs.boundaries.list, bname)
    end



    return nothing
end


