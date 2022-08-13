#----Definitions-----------------------

abstract type Unit end

struct UnitOp <: Unit
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

    # Number of historic data points
    numdata::Integer

    streamlist::StreamHistoryList
    # Connected streams with history
    inlets::Array{String, 1}
    outlets::Array{String, 1}

    # Inner constructor to validate input data
    function UnitOpHistory(name, numdata, streamlist, inlets, outlets)
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

        new(name, numdata, streamlist, inlets, outlets)
    end
end


struct UnitOpList
    list
end


struct UnitOpHistoryList
    list
end


#----Constructors----------------------


function UnitOp(name::String, streamlist::StreamList, inlets::Vector{String}, outlets::Vector{String})
    return UnitOp(name, streamlist, inlets, outlets, passive, nothing)
end


function UnitOp(name::String, streamlist::StreamList, inlets::Vector{String}, outlets::Vector{String}, f!::Function)
    return UnitOp(name, streamlist, inlets, outlets, f!, nothing)
end


"""
    UnitOpHistory(name, streamlist, inlets, outlets)

Constructor for a UnitOpHistory that defines the stream name and connected StreamHistory objects.
The internal field *numdata* is calculated from the connected streams.
"""
function UnitOpHistory(name, streamlist, inlets, outlets)
    numdata = length(streamlist[inlets[1]].totalmassflow)

    UnitOpHistory(name, numdata, streamlist, inlets, outlets)
end


function UnitOpList()
    l = Dict{String, UnitOp}()
    return UnitOpList(l)
end


function UnitOpHistoryList()
    l = Dict{String, UnitOpHistory}()
    return UnitOpHistoryList(l)
end


#----Active UnitOps--------------------


function passive(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)
    return nothing
end


function mixer!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)
    @assert length(outlets) == 1 "mixers can have only one outlet stream"

    old = streamlist[outlets[1]]
    oldname = old.name
    
    tempstream = Stream(old.name, old.complist, String[], Float64[])
    for inlet in inlets[1:end]
        tempstream += streamlist[inlet]
    end

    streamlist[oldname] = tempstream

    return nothing
end


#----Base overloads--------------------


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
    println(io, "Unit Operation: $(u.name)\n")
    println(io, "Feed streams:    ", [s for s in u.inlets])
    println(io)
    println(io, "Product streams: ", [s for s in u.outlets])
end


# Pretty printing for UnitOp objects
function Base.show(io::IO, u::UnitOpHistory)
    println(io, "Unit Operation:  $(u.name)\n")
    println(io, "Feed streams:    ", [s for s in u.inlets])
    println(io)
    println(io, "Product streams: ", [s for s in u.outlets])
    println(io)
    println(io, "Data length:     $(u.numdata)")
end