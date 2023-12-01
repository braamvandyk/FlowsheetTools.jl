# TODO Add uncertainty / variance to each stream for various flows

#----------------------------------------------------------------------------
#
#----Definitions-------------------------------------------------------------
#
#----------------------------------------------------------------------------


struct Stream
    name::String

    # Number of historic data points
    numdata::Integer
    
    complist::ComponentList
    
    # Flow data, mass, molar and atoms
    massflows::TimeArray
    moleflows::TimeArray
    totalmassflow::TimeArray
    atomflows::TimeArray
end


struct StreamList
    list::OrderedDict{String, Stream}
end


#----------------------------------------------------------------------------
#
#----Constructors------------------------------------------------------------
#
#----------------------------------------------------------------------------


"""
    function Stream(name, complist, comps, timestamps, flowdata, ismoleflow=false)

Constructor for a stream history object that defines the stream name and component flows
for various past measurements.

The mass flows are specified and molar flows and atomic molar flows calculated
    OR
The mole flows are specified and mass flows and atomic molar flows calculated

The flow data is passed as a matrix where each column represents a component. Internally, 
flows for all the components in the specified component list is stored, but only a sublist
need be specified in the constructor. Other components are assigned zero flows. This is to
expidite the addition of streams using the internal TimeArray objects.

It is recommended to rather create StreamHistory objects via readstreamhistory().
"""
function Stream(name, complist, comps, timestamps, flowdata, ismoleflow=false)
    numcomps = length(comps)
    numdata = size(flowdata, 1)

    
    # Sanity checks on the data
    size(flowdata, 2) != numcomps && error("mismatch between number of components and available data.")
    length(timestamps) != numdata && error("length mismatch between data and timestamps.")
    
    # Build the TimeArrays
    allcomps = names(complist) # Also those not present in the stream
    massflows = zeros(numdata, length(allcomps))
    moleflows = zeros(numdata, length(allcomps))

    # Collect the atoms to create a blank dictionary - all values are zero
    # This Dict only carries atoms actually present in the stream
    # When adding streams, this list is recreated from the new list of components
    emptyatomflows = Dict{String, Float64}()
    for compname in allcomps
        comp = complist[compname]
        for atom in comp.atoms
            if atom ∉ keys(emptyatomflows)
                emptyatomflows[atom] = 0.0
            end
        end
    end
    

    atomflows = Vector{Dict{String, Float64}}(undef, numdata)
    totalmassflows = zeros(numdata)
    
    for datum in 1:numdata
        # Get a copy of the empty atomic flows list
        _atomflows = deepcopy(emptyatomflows)

        # Here we loop through all the components in the complist and
        # use zeros when they are not present
        for (i, compname) in enumerate(allcomps)
            if compname in comps
                comp = complist[compname]
                idx = findfirst(x->x==compname, comps)

                if !ismoleflow
                    massflows[datum, i] = flowdata[datum, idx]
                    moleflows[datum, i] = flowdata[datum, idx]/comp.Mr
                else
                    moleflows[datum, i] = flowdata[datum, idx]
                    massflows[datum, i] = flowdata[datum, idx]*comp.Mr
                end

                totalmassflows[datum] += massflows[datum, i]
                
                for (j, atom) in enumerate(comp.atoms)
                    _atomflows[atom] += moleflows[datum, i]*comp.counts[j]
                end
            else
                massflows[datum, i] = 0.0
                moleflows[datum, i] = 0.0
            end
        end
        atomflows[datum] = _atomflows
    end

    massflows_ta = TimeArray(timestamps, massflows, allcomps)
    moleflows_ta = TimeArray(timestamps, moleflows, allcomps)
    totalmassflows_ta = TimeArray(timestamps, totalmassflows, [:totalmassflows])
    atomflows_ta = TimeArray(timestamps, atomflows, [:atomflows])
    

    return Stream(name, numdata, complist, massflows_ta, moleflows_ta, totalmassflows_ta, atomflows_ta)
end



"""
    Streamlist()

Constructor for an empty stream list. Streams are added when created via the @stream macro or
when read from file with readstreamhistory()
"""
function StreamList()
    l = OrderedDict{String, Stream}()
    return StreamList(l)
end


#----------------------------------------------------------------------------
# 
#----Base overloads----------------------------------------------------------
# 
#----------------------------------------------------------------------------


function Base.setindex!(A::StreamList, X::Stream, idx::String)
    if length(A.list) == 0
        A.list[idx] = X
    else
        # Verify that all the streams reference the same ComponentList.
        # Get the first item in the `list` field. As `list` is a `Dict`, this returns a `Pair`.
        # We get the value entry using the `second` field of the `Pair`, which returns a `Stream`,
        # of which we get the `complist` field.
        currentlist = first(A.list).second.complist
        current_ts = timestamp(first(A.list).second.massflows)

        X.complist != currentlist && error("all streams in StreamList must reference the same ComponentList")
        !all(timestamp(X.massflows) .== current_ts) && error("all streams in StreamList must have the same timestamps")

        A.list[idx] = X
    end

    return nothing
end


function Base.getindex(A::StreamList, idx::String)
    return A.list[idx]
end


function Base.getindex(A::StreamList, idxs::Vector{String})
    res = Stream[]
    for idx in idxs
        push!(res, A.list[idx])
    end
    return res
end


function Base.length(A::StreamList)
    return length(A.list)
end


"""
    Base.:+(a::Stream, b::Stream)

Extend the addition operator to add to streams to each other - a mixer.
All streams must refer to the same component list .
"""
function Base.:+(a::Stream, b::Stream)
    # Make sure the streams use the same system components and timestamps 
    # TimeArrays will add only matching timestamps, but this result in varying data lengths, which must be avoided
    a.complist != b.complist && error("cannot add streams with different system component lists")
    a.numdata != b.numdata && error("cannot add streams with different data lengths")
    !all(timestamp(a.massflows) .== timestamp(b.massflows)) && error("cannot add streams with different timestamps")

    comps = string.(colnames(a.massflows))
    timestamps = timestamp(a.massflows)
    flowdata = values(a.massflows .+ b.massflows)

    
    return Stream(a.name * "+" * b.name, a.complist, comps, timestamps, flowdata)
end


"""
    Base.*(a::T, b::Stream) where T <: Real

Extend the multiplication operator to scale a stream's flows by a scalar value.
Used in mass balance reconciliations to apply flow corrections.
"""
function Base.:*(a::T, b::Stream) where T <: Real
    comps = string.(colnames(b.massflows))
    timestamps = timestamp(b.massflows)
    flowdata = values(a .* b.massflows)
   
    return Stream(b.name, b.complist, comps, timestamps, flowdata)
end


function Base.:*(b::Stream, a::T) where T <: Real
    comps = string.(colnames(b.massflows))
    timestamps = timestamp(b.massflows)
    flowdata = values(a .* b.massflows)
   
    return Stream(b.name, b.complist, comps, timestamps, flowdata)
end


function  Base.copy(A::Stream)
    return deepcopy(A)
end


# # Pretty printing for stream objects
function Base.show(io::IO, stream::Stream)
    println(io, "Stream: $(stream.name)")
    println(io)

    if stream.numdata == 1
        header = string.(colnames(stream.massflows))
        pushfirst!(header, " ")
        massflows = ["Mass flows " (values(stream.massflows))]
        moleflows = ["Molar flows" (values(stream.moleflows))]
        data = vcat(massflows, moleflows)
        pretty_table(io, data, header = header)
        println(io)
        println(io, "Total mass flow: $(prettyround(values(stream.totalmassflow)[begin]))")
        println(io)
        pretty_table(io, values(stream.atomflows)[1], header = ["Atom", "Molar flow"])
    else
        println(io, "Mass flows:")
        pretty_table(io, stream.massflows, display_size=(14, -1))
        println(io)
        println(io, "Molar flows:")
        pretty_table(io, stream.moleflows, display_size=(14, -1))
        println(io)
        println(io)
        println(io, "Atom Flows:")
        pretty_table(io, stream.atomflows, display_size=(14, -1))
        println(io)
        println(io, "Data length:\t$(stream.numdata)")
        println(io, "Data starts:\t$(timestamp(stream.massflows)[begin])")
        println(io, "Data ends:\t$(timestamp(stream.massflows)[end])")   
    end
end


function  Base.show(io::IO, sl::StreamList)
    println(io, "Stream list:")
    println(io)
    println(io, "Streams:")
    if length(sl.list) > 0
        for (name, _) in sl
            println(io, "  ", name)
        end

        stream = first(sl.list).second
        
        println(io)
        println(io, "Components:")
        for (name, _) in stream.complist
            println(io, "  ", name)
        end

        println(io)
        println(io, "Data length:\t$(stream.numdata)")
        println(io, "Data starts:\t$(timestamp(stream.massflows)[begin])")
        println(io, "Data ends:\t$(timestamp(stream.massflows)[end])")    
    else
        println(io, "\tEmpty list")
    end
end


function Base.iterate(A::StreamList)
    return iterate(A.list)
end


function Base.iterate(A::StreamList, state)
    return iterate(A.list, state)
end


# ----------------------------------------------------------------------------
# 
#----Macros-------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------


"""
Defines a Stream() with the specified name and component mass flows and
add it to `sysstreams::StreamList`.

The names of components must match those in the specified componentlist
(`syscomps::ComponentList` in the example). 


    @stream mass begin
        "Ethylene" --> 2.8053
        "Ethane" --> 27.06192
        "Hydrogen" --> 2.21738
    end "Test" syscomps sysstreams

    @stream mole begin 
        "Ethylene" --> 0.1
        "Ethane" --> 0.9
        "Hydrogen" --> 1.1
    end "Product" syscomps sysstreams 

The created Stream object will have data at a single timestamp of DateTime(0).
"""
macro stream(flowtype::Symbol, ex::Expr, name::String, complist::Symbol, streamlist::Symbol)      
    local comps = String[]
    local flows = Float64[]
    local timestamps = [DateTime(0)]

    if flowtype == :mass
        local ismoleflow = false
    elseif flowtype == :mole
        local ismoleflow = true
    else
        error("flow basis specification must be \"mass\" or \"mole\"")
    end

    for line in ex.args
        match_comp = @capture(line, comp_ --> flow_)
        if match_comp
            local component = eval(comp)
            local flowval = eval(flow)
            if any(x -> x == component, comps)
                i = findfirst(x -> x == component, comps)
                flows[i] += flowval
            else
                push!(comps, component)
                push!(flows, flowval)
            end
        end
    end
    return :($(esc(streamlist))[$name] = Stream($name, $(esc(complist)), $comps, $timestamps, $(flows'), $ismoleflow))
end


#----------------------------------------------------------------------------
# 
#----Utilities---------------------------------------------------------------
# 
#----------------------------------------------------------------------------


"""
    function copystream!(list::StreamList, from::String, to::String)

Copy a stream in the stream list
"""
function copystream!(list::StreamList, from::String, to::String; factor=1.0)   
    fromstream = list[from]
    comps = string.(colnames(fromstream.massflows))
    timestamps = timestamp(fromstream.massflows)
    flowdata = factor .* values(fromstream.massflows)
   
    list[to] = Stream(to, fromstream.complist, comps, timestamps, flowdata)

    return nothing
end


"""
    function deletestream!(list::StreamList, from::String)

Delete a stream from the stream list
"""
function deletestream!(list::StreamList, from::String)
    delete!(list.list, from)

    return nothing
end



"""
    function renamestream!(list::StreamList, from::String, to::String)

Rename a stream in the stream list
"""
function renamestream!(list::StreamList, from::String, to::String)
    fromstream = list[from]
    comps = string.(colnames(fromstream.massflows))
    timestamps = timestamp(fromstream.massflows)
    flowdata = values(fromstream.massflows)

    list[to] = Stream(to, fromstream.complist, comps, timestamps, flowdata)
    deletestream!(list, from)
    return nothing
end

"""
    renamestream(strm::Stream, to::String)

Return a copy of the stream with a new name
"""
function renamestream(fromstream::Stream, to::String)
    comps = string.(colnames(fromstream.massflows))
    timestamps = timestamp(fromstream.massflows)
    flowdata = values(fromstream.massflows)

    return Stream(to, fromstream.complist, comps, timestamps, flowdata)
end


"""
    emptystream(list::StreamList, name::String)

Returns a stream with the same components and timestamp as the specified StreamList and all flows set to zero. 
"""
function emptystream(list::StreamList, name::String)
    # Take the first stream in the StreamList as a reference
    refstrm = first(list).second

    return Stream(name, refstrm.complist, string.(colnames(refstrm.massflows)), timestamp(refstrm.massflows), zeros(size(refstrm.massflows)))
end


"""

    readstreamhistory(filename, streamname, complist)

Reads in a stream history file (CSV file).

The data should be in the format of TimeStamp (yyyy/mm/dd HH:MM) in first column,
then each subsequent column holding a component's mass flows, with the heading the
name of the component.

**The names should match those in the system component list that is passed to the function.**

"""
function readstreamhistory(filename, streamname, complist; ismoleflow=false)

    data, header = readdlm(filename, ',', '\n', header=true)
    comps = string.(header[2:end]) # readdlm returns an array of AbstractStrings for some reason
    flows = data[:, 2:end]
    timestamps = DateTime.(data[:, 1], "yyyy/mm/dd HH:MM")

    return Stream(streamname, complist, comps, timestamps, flows, ismoleflow)
end


"""
    refreshcomplist(streamlist::StreamList)

Refresh all the streams in the StreamList to reflect components added to the StreamList after the streams were created.
"""
function refreshcomplist(streamlist::StreamList)
    for (name, stream) in streamlist
        complist = stream.complist
        comps = string.(colnames(stream.massflows))
        timestamps = timestamp(stream.massflows)
        flowdata = values(stream.massflows)

        streamlist[name] = Stream(name, complist, comps, timestamps, flowdata)
    end
end
