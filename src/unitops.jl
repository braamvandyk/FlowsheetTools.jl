#----------------------------------------------------------------------------
#
#----Definitions-------------------------------------------------------------
#
#----------------------------------------------------------------------------


struct UnitOp
    name::String

    streamlist::StreamList

    # Connected streams
    inlets::Array{String, 1}
    outlets::Array{String, 1}

    f!::Function
    params
end

function (u::UnitOp)(newparams = nothing)
    if isnothing(newparams)
        u.f!(u.streamlist, u.outlets, u.inlets, u.params)
    else
        u.f!(u.streamlist, u.outlets, u.inlets, newparams)
    end
end


struct UnitOpHistory
    name::String

    streamlist::StreamHistoryList
    # Connected streams with history
    inlets::Array{String, 1}
    outlets::Array{String, 1}

    f!::Function
    params

    # Inner constructor to validate input data
    function UnitOpHistory(name, streamlist, inlets, outlets, f!, params)
        numdata = length(streamlist[inlets[1]].totalmassflow)
        timestamps = streamlist[inlets[1]].timestamps
        
        for inletname in inlets
            inlet = streamlist[inletname]

            inlet.numdata != numdata && error("all streams must must have similar history lengths.")
            for j in eachindex(timestamps)
                timestamps[j] != inlet.timestamps[j] && error("all stream histories must have matching timestamps.")
            end
        end
        for outletname in outlets
            outlet = streamlist[outletname]

            outlet.numdata != numdata && error("all streams must must have similar history lengths.")
            for j in eachindex(timestamps)
                timestamps[j] != outlet.timestamps[j] && error("all stream histories must have matching timestamps.")
            end
        end

        new(name, streamlist, inlets, outlets, f!, params)
    end
end


function (u::UnitOpHistory)(newparams = nothing)
    if isnothing(newparams)
        u.f!(u.streamlist, u.outlets, u.inlets, u.params)
    else
        u.f!(u.streamlist, u.outlets, u.inlets, newparams)
    end
end


struct UnitOpList
    list
end


struct UnitOpHistoryList
    list
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
    UnitOpHistory(name::String, streamlist::StreamHistoryList, inlets::Vector{String}, outlets::Vector{String})
    UnitOpHistory(name::String, streamlist::StreamHistoryList, inlets::Vector{String}, outlets::Vector{String}, f!::Function)
    UnitOpHistory(name::String, streamlist::StreamHistoryList, inlets::Vector{String}, outlets::Vector{String}, f!::Function, params)

Constructor for a UnitOp. It is recommended to rather use the @unitop macro.
"""
function UnitOpHistory(name::String, streamlist::StreamHistoryList, inlets::Vector{String}, outlets::Vector{String})
    return UnitOpHistory(name, streamlist, inlets, outlets, passive, nothing)
end


function UnitOpHistory(name::String, streamlist::StreamHistoryList, inlets::Vector{String}, outlets::Vector{String}, f!::Function)
    return UnitOpHistory(name, streamlist, inlets, outlets, f!, nothing)
end


"""
    UnitOpList()

Constructor for an empty UnitOpList. Unit operations are added when created with the @unitop macro.
"""
function UnitOpList()
    l = Dict{String, UnitOp}()
    return UnitOpList(l)
end


"""
    UnitOpHistoryList()

Constructor for an empty UnitOpHistoryList. Unit operations are added when created with the @unitophist macro.
"""
function UnitOpHistoryList()
    l = Dict{String, UnitOpHistory}()
    return UnitOpHistoryList(l)
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
function passive(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)
    return nothing
end


"""
    mixer!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)
    mixer!(streamlist::StreamHistoryList, outlets::Vector{String}, inlets::Vector{String}, params)

Calculation for mixer UnitOps. Combines all feed streams into a single outlet stream.
Will error if more than one outlet stream is defined. Values in `params` are ignored.
"""
function mixer!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)
    @assert length(outlets) == 1 "mixers can have only one outlet stream"

    old = streamlist[outlets[1]]
        
    tempstream = Stream(old.name, old.complist, String[], Float64[])
    for inlet in inlets[1:end]
        tempstream += streamlist[inlet]
    end

    tempstream = renamestream(tempstream, old.name)
    streamlist[old.name] = tempstream
    
    return nothing
end


function mixer!(streamlist::StreamHistoryList, outlets::Vector{String}, inlets::Vector{String}, params)
    @assert length(outlets) == 1 "mixers can have only one outlet stream"

    old = streamlist[inlets[1]]
        
    tempstream = StreamHistory(old.name, old.complist, String[] DateTime[], Matrix{Float64}())
    for inlet in inlets[1:end]
        tempstream += streamlist[inlet]
    end

    tempstream = renamestream(tempstream, old.name)
    streamlist[old.name] = tempstream
    
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


function Base.setindex!(A::UnitOpHistoryList, X::UnitOpHistory, idx::String)
    if length(A.list) == 0
        A.list[idx] = X
    else
        # Verify that all the streams reference the same ComponentList.
        # Get the first item in the `list` field. As `list` is a `Dict`, this returns a `Pair`.
        # We get the value entry using the `second` field of the `Pair`, which returns a `UnitOpHistory`,
        # of which we get the `streamlist` field.
        currentlist = first(A.list).second.streamlist
        X.streamlist != currentlist && error("all streams in UnitOpHistoryList must reference the same StreamHistoryList")

        A.list[idx] = X
    end
end


function Base.getindex(A::UnitOpHistoryList, idx::String)
    return A.list[idx]
end


# Pretty printing for UnitOp objects
function Base.show(io::IO, u::UnitOp)
    println(io, "Unit Operation:  $(u.name)")
    println(io, "Feed streams:    ", [s for s in u.inlets])
    println(io)
    println(io, "Product streams: ", [s for s in u.outlets])
end


# Pretty printing for UnitOp objects
function Base.show(io::IO, u::UnitOpHistory)
    println(io, "Unit Operation:  $(u.name)")
    println(io, "Feed streams:    ", [s for s in u.inlets])
    println(io)
    println(io, "Product streams: ", [s for s in u.outlets])
    println(io)
    strm = first(u.streamlist.list)
    println(io, "Data length:\t$(strm.second.numdata)")
    println(io, "Data starts:\t$(strm.second.timestamps[begin])")
    println(io, "Data ends:\t$(strm.second.timestamps[end])") 
end


function  Base.show(io::IO, uol::UnitOpList)
    println(io, "UnitOp list:")
    for unitop in uol.list
        println(io, "  ", unitop.first)
    end
end


function  Base.show(io::IO, uol::UnitOpHistoryList)
    println(io, "UnitOpHistory list:")
    if length(uol.list) > 0
        for unitop in uol.list
            println(io, "  ", unitop.first)
        end
        println(io)
        sl = first(uol.list).second.streamlist
        inlets = first(uol.list).second.inlets
        strm = sl[inlets[1]]
        println(io, "Data length:\t$(strm.numdata)")
        println(io, "Data starts:\t$(strm.timestamps[begin])")
        println(io, "Data ends:\t$(strm.timestamps[end])")         
    else
        println(io, "Empty list")
    end
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


"""
Defines a `UnitOpHist with` the specified name and connected streams.

@unitophist begin
    inlets --> ["H2", "C2"]
    outlets --> ["Mixed"]
    calc --> mixer!
end "Mixer" histstreams histunitops

This will create a UnitOpHist named "Mixer", saved into histunitops["Mixer"].
The specified streams refer to entries in histstreams::StreamHistoryList.

Two optional parameters may be specified, i.e. calc and params:

calc is a function, e.g. 
    function mixer!(streamlist::StreamHistoryList, outlets::Vector{String}, inlets::Vector{String}, params)
params is an iterable passed to the function.

Together, these specify a calculation to perform to calculate the outlet(s) from the inlet(s), when
calling `histunitops["Mixer"]()`

Alternative parameters may also be specified in the call, e.g.
    histunitops["Mixer"]([1, 2, 3])

"""
macro unitophist(ex::Expr, name::String, streamlist::Symbol, unitoplist::Symbol)      
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
        return :($(esc(unitoplist))[$name] = UnitOpHistory($name, $(esc(streamlist)), $inlets, $outlets))
    else    
        return :($(esc(unitoplist))[$name] = UnitOpHistory($name, $(esc(streamlist)), $inlets, $outlets, $(esc(func)), $(esc(params))))
    end
end