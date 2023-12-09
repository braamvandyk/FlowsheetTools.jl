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

            inlet.numdata != numdata && error("all streams must must have identical history lengths.")
            !all(inlet_ts .== timestamps) && error("all streams must must have identical timestamps.")
        end
        for outletname in outlets
            outlet = streamlist[outletname]
            outlet_ts = timestamp(outlet.massflows)

            outlet.numdata != numdata && error("all streams must must have identical history lengths.")
            !all(outlet_ts .== timestamps) && error("all streams must must have identical timestamps.")
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


"""
    mixer!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)

Calculation for mixer UnitOps. Combines all feed streams into a single outlet stream.
Will error if more than one outlet stream is defined. Values in `params` are ignored.
"""
function mixer!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)
    @assert length(outlets) == 1 "mixers can have only one outlet stream"

    name = streamlist[outlets[1]].name     
    tempstream = sum(streamlist[inlets])

    tempstream = renamestream(tempstream, name)
    streamlist[name] = tempstream
    
    return nothing
end


"""
    flowsplitter!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)

Calculation for flowsplitter UnitOps. Splits combined feed stream according to params - an iterable with 
split fractions to (n - 1) streams. The last stream gets the balance.

"""
function flowsplitter!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)
    numouts = length(outlets)
    sumfracs = sum(params)

    @assert length(params) == (numouts - 1) "incorrect number of fractions - must be one less than outlets"
    @assert 0.0 <= sumfracs <= 1.0 "invalid split fractions specified"

    totalin = sum(streamlist[inlets])
    
    for (i, outname) in enumerate(outlets)
        if i < numouts
            frac = params[i]
        else
            frac = 1.0 - sumfracs
        end

        streamlist[outname] = renamestream(frac * totalin, outname)
    end
    
    return nothing
end


"""
    
    componentplitter!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)

Calculation for componentplitter UnitOps. Splits combined feed stream according to params - a nested Dict with 
split fractions for each component to (n - 1) streams. The last stream gets the balance.

Example:
    params = Dict([
    "Hydrogen" => Dict(["Product1" => 0.5]),
    "Ethylene" => Dict(["Product2" => 0.3])
    ])

"""
function componentplitter!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)
    complist = first(streamlist).second.complist
    numcomps = length(complist)
    numouts = length(outlets)
    totalin = sum(streamlist[inlets])
    fractions = fill(-1.0, (numcomps, numouts))
    compnames = Array{String}(undef, numcomps)

    # Loop through the specified component splits and assign the specified fractions to their streams
    # Any fraction not specified will remain as -1.0, which will be processed in the next step
    for (compindex, comp) in enumerate(complist)
        compname = comp.first
        compnames[compindex] = compname
        if compname in keys(params)
            #  This component has specified splits
            splits = params[compname]
            @assert length(splits) == numouts - 1 "incorrect number of fractions specified for component $compname"
            for split in splits
                stream = split.first
                @assert stream in outlets "specified stream not in list of outlet streams"
                streamindex = findfirst(str -> str == stream, outlets)
                fraction = split.second
                @assert 0.0 <= fraction <= 1.0 "invalid fraction specified for $compname in stream @stream"
                fractions[compindex, streamindex] = fraction
            end
        end
    end
    
    #  Now we look for unspecified fractions and set those.
    for row in 1:size(fractions)[1]
        sumfracs = sum(fractions[row, :])
        if sumfracs == -numouts
            # No splits specified for this component -> everything goes to last stream
            fractions[row, :] .= 0.0
            fractions[row, end] = 1.0
        else
            # We already know from the previous step that n - 1 streams are assigned for
            # each component with specified splits, or the assert would have failed.
            # This means the only unspecified component still has a value of -1.0. This
            # fraction must be the 1 - sum(specified fractions)
            for j in 1:numouts
                if fractions[row, j] == -1.0
                    fractions[row, j] = 1.0 - (sumfracs + 1.0)
                    break
                end
            end
        end

    end

    total_massflows = values(totalin.massflows)
    timestamps = timestamp(totalin.massflows)
    for streamindex in 1:numouts
        massflows = copy(total_massflows)
        for datum in 1:size(total_massflows)[1]
            massflows[datum, :] .*= fractions[:, streamindex]
        end
        streamlist[outlets[streamindex]] = Stream(outlets[streamindex], complist, compnames, timestamps, massflows)
    end

    return nothing
end


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
        X.streamlist != currentlist && error("all unit operations in UnitOpList must reference the same StreamList")

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
