# TODO Add a splitter unit
# TODO Add a stoichiometric reactor with multiple reactions and conversions

#----------------------------------------------------------------------------
#
#----Definitions-------------------------------------------------------------
#
#----------------------------------------------------------------------------


struct UnitOp
    name::String

    streamlist::StreamList
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
        u.f!(u.streamlist, u.outlets, u.inlets, u.params)
    else
        u.f!(u.streamlist, u.outlets, u.inlets, newparams)
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
        currentlist = first(A.list).second.streamlist
        X.streamlist != currentlist && throw(ArgumentError("all unit operations in UnitOpList must reference the same StreamList"))

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
    
    strm = first(u.streamlist).second
    
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
        sl = uo.streamlist
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
end "Mixer" sysstreams sysunitops

This will create a UnitOp named "Mixer", saved into sysunitops["Mixer"].
The specified streams refer to entries in sysstreams::StreamList.

Two optional parameters may be specified, i.e. calc and params:

calc is a function, e.g. 

    function mixer!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)

params is an iterable passed to the function.

Together, these specify a calculation to perform to calculate the outlet(s) from the inlet(s), when
calling `sysunitops["Mixer"]()`

Alternative parameters may also be specified in the call, e.g.
    sysunitops["Mixer"]([1, 2, 3])

"""
macro unitop(ex::Expr, name::String, streamlist::Symbol, unitoplist::Symbol)      
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
        return :($(esc(unitoplist))[$name] = UnitOp($name, $(esc(streamlist)), $inlets, $outlets))
    else    
        return :($(esc(unitoplist))[$name] = UnitOp($name, $(esc(streamlist)), $inlets, $outlets, $(esc(func)), $(esc(params))))
    end
end
