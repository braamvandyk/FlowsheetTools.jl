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
    
    comps::ComponentList
    
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
    size(flowdata, 2) != numcomps && throw(DimensionMismatch("mismatch between number of components and available data."))
    length(timestamps) != numdata && throw(DimensionMismatch("length mismatch between data and timestamps."))
    
    # Build the TimeArrays
    allcomps = members(complist) # Also those not present in the stream
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

    StreamList()

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
        currentlist = first(A.list).second.comps
        current_ts = timestamp(first(A.list).second.massflows)

        X.comps != currentlist && throw(ArgumentError("all streams in StreamList must reference the same ComponentList"))
        !all(timestamp(X.massflows) .== current_ts) && throw(ArgumentError("all streams in StreamList must have the same timestamps"))

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


function Base.length(A::Stream)
    return A.numdata
end


function Base.:+(a::Stream, b::Stream)
    # Extend the addition operator to add to streams to each other - a mixer.
    # All streams must refer to the same component list.

    # Make sure the streams use the same system components and timestamps 
    # TimeArrays will add only matching timestamps, but this result in varying data lengths, which must be avoided
    a.comps != b.comps && throw(ArgumentError("cannot add streams with different system component lists"))
    a.numdata != b.numdata && throw(DimensionMismatch("cannot add streams with different data lengths"))
    !all(timestamp(a.massflows) .== timestamp(b.massflows)) && throw(ArgumentError("cannot add streams with different timestamps"))

    comps = string.(colnames(a.massflows))
    timestamps = timestamp(a.massflows)
    flowdata = values(a.massflows .+ b.massflows)

    
    return Stream(a.name * "+" * b.name, a.comps, comps, timestamps, flowdata)
end


# Extend the multiplication operator to scale a stream's flows by a scalar value.
# Used in mass balance reconciliations to apply flow corrections.
function Base.:*(a::T, b::Stream) where T <: Real
    comps = string.(colnames(b.massflows))
    timestamps = timestamp(b.massflows)
    flowdata = values(a .* b.massflows)
   
    return Stream(b.name, b.comps, comps, timestamps, flowdata)
end


function Base.:*(b::Stream, a::T) where T <: Real
    comps = string.(colnames(b.massflows))
    timestamps = timestamp(b.massflows)
    flowdata = values(a .* b.massflows)
   
    return Stream(b.name, b.comps, comps, timestamps, flowdata)
end


function Base.:≈(a::Stream, b::Stream)
    # Extend the ≈ operator to check if two streams have approximately equal flows.
    # The check is done internally on molar flows.
    return all(values(a.moleflows) .≈ values(b.moleflows))
end


function Base.:(==)(a::Stream, b::Stream)
    # Extend the == operator to check if two streams have equal flows.
    # The check is done internally on molar flows.
    # Since flows are floating point values, it is recommended to rather use ≈ for comparisons.
    # Will also insist on identical names!
    return (all(values(a.moleflows) .== values(b.moleflows))) && (a.name == b.name)
end


function  Base.copy(A::Stream)
    return deepcopy(A)
end


function Base.show(io::IO, stream::Stream)
    # Pretty printing for stream objects
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
        for (name, _) in stream.comps
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
add it to `fs.streams::StreamList`.

The names of components must match those in the specified componentlist
(`fs.comps::ComponentList` in the example). 


    @stream mass begin
        "Ethylene" --> 2.8053
        "Ethane" --> 27.06192
        "Hydrogen" --> 2.21738
    end "Test" fs

    @stream mole begin 
        "Ethylene" --> 0.1
        "Ethane" --> 0.9
        "Hydrogen" --> 1.1
    end "Product" fs 

The created Stream object will have data at a single timestamp of DateTime(0).

"""
macro stream(flowtype::Symbol, ex::Expr, name::String, fs::Symbol)      
    local comps = String[]
    local flows = Float64[]
    local timestamps = [DateTime(0)]

    if flowtype == :mass
        local ismoleflow = false
    elseif flowtype == :mole
        local ismoleflow = true
    else
        throw(ArgumentError("flow basis specification must be \"mass\" or \"mole\""))
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
    return :($(esc(fs)).streams[$name] = Stream($name, $(esc(fs)).comps, $comps, $timestamps, $(flows'), $ismoleflow))
end


#----------------------------------------------------------------------------
# 
#----Utilities---------------------------------------------------------------
# 
#----------------------------------------------------------------------------


"""

    function copystream!(fs, from, to; factor=1.0)

Copy a stream in the stream list of fs::Flowsheet
`from` and `to` are the names of the source and destination streams
The source stream's flows are multiplied by `factor` before being copied

"""
function copystream!(fs, from, to; factor=1.0)
    @argcheck fs isa Flowsheet "fs must be a Flowsheet"

    fromstream = fs.streams[from]
    comps = string.(colnames(fromstream.massflows))
    timestamps = timestamp(fromstream.massflows)
    flowdata = factor .* values(fromstream.massflows)
   
    fs.streams[to] = Stream(to, fs.comps, comps, timestamps, flowdata)

    return nothing
end


"""

    function deletestream!(fs, name)

Delete a stream from the Flowsheet's StreamList.
"""
function deletestream!(fs, name)
    @argcheck fs isa Flowsheet "fs must be a Flowsheet"
    if name in keys(fs.streams.list)
        delete!(fs.streams.list, name)

        # Also clear the unitops to ensure consistency
        for (uname, unitop) in fs.unitops.list
            if name in unitop.inlets || name in unitop.outlets
                deleteunitop!(fs, uname) # Will cascase delete boundaries as well
            end
        end
    end

    return nothing
end


"""

function deletestreams!(fs)

Delete all streams from the Flowsheet's StreamList.

"""
function deletestreams!(fs)
    @argcheck fs isa Flowsheet "fs must be a Flowsheet"

    for name in keys(fs.streams.list)
        delete!(fs.streams.list, name)
    end

    # Also kill all the streams and boundaries
    deleteunitops!(fs)

    return nothing
end



"""

    function renamestream!(fs, from, to)

Rename a stream in the stream list in fs::Flowsheet
`from` and `to` are the old and new names of the stream

"""
function renamestream!(fs, from, to)
    fromstream = fs.streams[from]
    comps = string.(colnames(fromstream.massflows))
    timestamps = timestamp(fromstream.massflows)
    flowdata = values(fromstream.massflows)

    fs.streams[to] = Stream(to, fromstream.comps, comps, timestamps, flowdata)
    deletestream!(fs, from)
    return nothing
end


"""

    renamestream(stream, newname)

Returns a new Stream object with the same information as the specificied stream, but a new name.

"""
function renamestream(stream, newname)
    comps = members(stream.comps)
    timestamps = timestamp(stream.massflows)
    flowdata = values(stream.massflows)

    return Stream(newname, stream.comps, comps, timestamps, flowdata)
end    

"""

    addemptystream!(fs, name)

Returns a stream with the same components and timestamp as the other streams in `fs::Flowsheet` with all flows set to zero. 

"""
function addemptystream!(fs, name::String)
    # Take the first stream in the StreamList as a reference to get the components and timestamps.
    # Here `first()` returns a `Pair` of a name and a stream, so we call `second()` to get the `Stream` object
    refstrm = try
        first(fs.streams).second
    catch
        throw(ArgumentError("fs must contain at least one stream"))
    end

    fs.streams[name] = Stream(name, fs.comps, string.(colnames(refstrm.massflows)), timestamp(refstrm.massflows), zeros(size(refstrm.massflows)))

    return nothing
end


"""

    addfixedstream!(fs, name, flows; ismoleflow=false)
    
Returns a stream with the same components and timestamp as the other streams in `fs::Flowsheet` with all flows set to a specified constant values.
Use the second form if no other streams exist, to specifiy the names of the components.

`flows` should contain the constant flowrates for all the components in the stream.

"""
function addfixedstream!(fs, name, flows; ismoleflow=false)
    # Take the first stream in the StreamList as a reference to get the components and timestamps.
    # Here `first()` returns a `Pair` of a name and a stream, so we call `second()` to get the `Stream` object
    refstrm = try
        first(fs.streams).second
    catch
        throw(ArgumentError("fs must contain at least one stream"))
    end

    flowdata = repeat(flows', outer = length(timestamp(refstrm.massflows)))
    fs.streams[name] = Stream(name, refstrm.comps, string.(colnames(refstrm.massflows)), timestamp(refstrm.massflows), flowdata, ismoleflow)

    return nothing
end


"""

     readstreamhistory!(fs, streamname, filename; ismoleflow=false)

Reads in a stream history file (CSV file) into the stream in fs.streams.

The data should be in the format of TimeStamp (yyyy/mm/dd HH:MM) in first column,
then each subsequent column holding a component's mass flows, with the heading the
name of the component.

**The names should match those in the specified component list that is passed to the function.**

"""
function readstreamhistory!(fs, streamname, filename; ismoleflow=false)

    data, header = readdlm(filename, ',', '\n', header=true)
    comps = string.(header[2:end]) # readdlm returns an array of AbstractStrings for some reason
    flows = data[:, 2:end]
    timestamps = DateTime.(data[:, 1], "yyyy/mm/dd HH:MM")

    fs.streams[streamname] = Stream(streamname, fs.comps, comps, timestamps, flows, ismoleflow)

    return nothing
end


"""

    writestreamhistory(filename, streamname)

Writes in a stream history file (CSV file).


"""
function writestreamhistory(stream, filename, moleflow=false)
    if moleflow
        writetimearray(stream.moleflows, filename)
    else
        writetimearray(stream.massflows, filename)
    end

    return nothing
end


"""

    writestreamhistories(fs, filenames, moleflow=false)

Writes stream history files (CSV files) for each stream in a Flowsheet, fs.
`filenames` is a dictionary of stream names and filenames, e.g. `filenames["stream1"] = "stream1.csv"`.

"""
function writestreamhistories(fs, filenames, moleflow=false)

    for (name, stream) in fs.streams
        writestreamhistory(stream, filenames[name], moleflow)
    end

    return nothing
end


"""
    refreshcomplist(fs)

Refresh all the streams in the StreamList to reflect components added to the StreamList after the streams were created.
"""
function refreshcomplist(fs)
    for (name, stream) in fs.streams
        complist = stream.comps
        comps = string.(colnames(stream.massflows))
        timestamps = timestamp(stream.massflows)
        flowdata = values(stream.massflows)

        fs.streams[name] = Stream(name, complist, comps, timestamps, flowdata)
    end
end

