#----------------------------------------------------------------------------
#
#----Definitions-------------------------------------------------------------
#
#----------------------------------------------------------------------------


struct Stream
    name::String

    complist::ComponentList
    comps::Array{String, 1}

    # Specifiy mass flows, calculate mole flows
    massflows::Array{Float64, 1}
    moleflows::Array{Float64, 1}
    totalmassflow::Float64

    # Molar flow of atoms, calculated from above
    atomflows::Dict{String, Float64}
end


struct StreamHistory
    name::String

    # Number of historic data points
    numdata::Integer
    
    complist::ComponentList
    comps::Array{String, 1}

    timestamps::Array{DateTime, 1}
    
    # Specifiy mass flows, calculate mole flows
    massflows::Array{Float64, 2}
    moleflows::Array{Float64, 2}
    totalmassflow::Array{Float64, 1}

    # Molar flow of atoms, calculated from above
    atomflows::Array{Dict{String, Float64}, 1}
end


struct StreamList
    list
end


struct StreamHistoryList
    list
end


#----------------------------------------------------------------------------
#
#----Constructors------------------------------------------------------------
#
#----------------------------------------------------------------------------


"""
    function Stream(name, complist, comps, flows, ismoleflow=false)

Constructor for a stream that defines the stream name and component flows.
The mass flows are specified and a molar composition and atomic molar flows calculated.
"""
function Stream(name, complist, comps, flows, ismoleflow=false)
    numcomps = length(comps)
    
    length(flows) != numcomps && error("mismatch between number of components and available data.")

    massflows = zeros(numcomps)
    moleflows = zeros(numcomps)

    atomflows = Dict{String, Float64}()
    totalmassflow = 0.0
    

    for (i, compname) in enumerate(comps)
        comp = complist[compname]

        if !ismoleflow
            massflows[i] = flows[i]
            moleflows[i] = massflows[i]/comp.Mr
        else
            moleflows[i] = flows[i]
            massflows[i] = moleflows[i]*comp.Mr
        end
            
        totalmassflow += massflows[i]
            
        for (j, atom) in enumerate(comp.atoms)
            if atom in keys(atomflows)
                atomflows[atom] += moleflows[i]*comp.counts[j]
            else
                atomflows[atom] = moleflows[i]*comp.counts[j]
            end
        end

    end

    return Stream(name, complist, comps, massflows, moleflows, totalmassflow, atomflows)
end


"""
    function StreamHistory(name, complist, comps, timestamps, flowshistory, ismoleflow=false)

Constructor for a stream history object that defines the stream name and component flows
for various past measurements.
The mass flows are specified and a molar composition and atomic molar flows calculated.
The mass flow history is passed as a matrix where each column is a datum
"""
function StreamHistory(name, complist, comps, timestamps, flowshistory, ismoleflow=false)
    numcomps = length(comps)
    numdata = size(flowshistory, 2)

    size(flowshistory, 1) != numcomps && error("mismatch between number of components and available data.")
    length(timestamps) != numdata && error("length mismatch between data and timestamps.")

    massflowshistory = zeros(numcomps, numdata)
    moleflowshistory = zeros(numcomps, numdata)

    atomflowshistory = Array{Dict{String, Float64}, 1}(undef, numdata)
    totalmassflowhistory = zeros(numdata)
    
    for datum in 1:numdata
        atomflows = Dict{String, Float64}()

        for (i, compname) in enumerate(comps)
            comp = complist[compname]

            if !ismoleflow
                massflowshistory[i, datum] = flowshistory[i, datum]
                moleflowshistory[i, datum] = flowshistory[i, datum]/comp.Mr
            else
                moleflowshistory[i, datum] = flowshistory[i, datum]
                massflowshistory[i, datum] = flowshistory[i, datum]*comp.Mr
            end

            totalmassflowhistory[datum] += massflowshistory[i, datum]
            
            for (j, atom) in enumerate(comp.atoms)
                if atom in keys(atomflows)
                    atomflows[atom] += moleflowshistory[i, datum]*comp.counts[j]
                else
                    atomflows[atom] = moleflowshistory[i, datum]*comp.counts[j]
                end
            end
        end
        atomflowshistory[datum] = atomflows
    end

    return StreamHistory(name, numdata, complist, comps, timestamps, massflowshistory, moleflowshistory,
                         totalmassflowhistory, atomflowshistory)
end


function StreamList()
    l = Dict{String, Stream}()
    return StreamList(l)
end


function StreamHistoryList()
    l = Dict{String, StreamHistory}()
    return StreamHistoryList(l)
end

#----------------------------------------------------------------------------
# 
#----Base overloads----------------------------------------------------------
# 
#----------------------------------------------------------------------------


"""
    Base.:+(a::Stream, b::Stream)

Extend the addition operator to add to streams to each other - a mixer.
It is assumed that the streams will have different components in arbitrary order.
"""
function Base.:+(a::Stream, b::Stream)
    # Make sure the streams use the same system components
    a.complist != b.complist && error("cannot add streams with different system component lists")
    
    comps = copy(a.comps)
    massflows = copy(a.massflows)

    compmap = indexin(b.comps, comps)
    for i in eachindex(b.comps)
        if isnothing(compmap[i])
            # Add the component
            push!(comps, b.comps[i])
            push!(massflows, b.massflows[i])
        else
            # Add the flow to the existing component
            massflows[compmap[i]] += b.massflows[i]
        end
    end

    return Stream(a.name * "+" * b.name, a.complist, comps, massflows)
end


"""
    Base.+(a::StreamHistory, b::StreamHistory)

Extend the addition operator to add to streams histories to each other - a mixer.
It is assumed that the streams histories will have different components in arbitrary order.
"""
function Base.:+(a::StreamHistory, b::StreamHistory)
    # Make sure the streams use the same system components
    a.complist != b.complist && error("cannot add streams with different system component lists")

    # Check that the data peridos are the same.
    # The constructor already checks that the timestamp length matches the data length.
    # We don't check the comps, as the streams could have different compositions, which we mix.
    a.numdata != b.numdata && error("history length of a not identical to that of b")

    # Check that the timestamps are identical!
    for i in eachindex(a.timestamps)
        a.timestamps[i] != b.timestamps[i] && error("timestamp values do not match at entry $i")
    end

    comps = copy(a.comps)
    massflows = copy(a.massflows)

    compmap = indexin(b.comps, comps)
    for i in eachindex(b.comps)
        if isnothing(compmap[i])
            # Add the component
            push!(comps, b.comps[i])
            vcat(massflows, b.massflows[i, :]')
        else
            # Add the flow to the existing component
            massflows[compmap[i], :] .+= b.massflows[i, :]
        end
    end

    return StreamHistory(a.name * "+" * b.name, a.complist, comps, a.timestamps, massflows)
end


"""
    Base.*(a::T, b::Stream) where T <: Real

Extend the multiplication operator to scale a stream's flows by a scalar value.
Used in mass balance reconciliations to apply flow corrections.
"""
function Base.:*(a::T, b::Stream) where T <: Real
    return Stream(b.name, b.complist, b.comps, a .* b.massflows)
end

function Base.:*(b::Stream, a::T) where T <: Real
    return Stream(b.name, b.complist, b.comps, a .* b.massflows)
end

"""
    Base.*(a::T, b::StreamHistory) where T <: Real

Extend the multiplication operator to scale a stream history's flows by a scalar value.
Used in mass balance reconciliations to apply flow corrections.
"""
function Base.:*(a::T, b::StreamHistory) where T <: Real
    return StreamHistory(b.name, b.complist, b.comps, b.timestamps, a .* b.massflows)
end

function Base.:*(b::StreamHistory, a::T) where T <: Real
    return StreamHistory(b.name, b.complist, b.comps, b.timestamps, a .* b.massflows)
end


function Base.setindex!(A::StreamList, X::Stream, idx::String)
    if length(A.list) == 0
        A.list[idx] = X
    else
        # Verify that all the streams reference the same ComponentList.
        # Get the first item in the `list` field. As `list` is a `Dict`, this returns a `Pair`.
        # We get the value entry using the `second` field of the `Pair`, which returns a `Stream`,
        # of which we get the `complist` field.
        currentlist = first(A.list).second.complist
        X.complist != currentlist && error("all streams in StreamList must reference the same ComponentList")

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


function Base.setindex!(A::StreamHistoryList, X::StreamHistory, idx::String)
    if length(A.list) == 0
        A.list[idx] = X
    else
        # Verify that all the streams reference the same ComponentList.
        # Get the first item in the `list` field. As `list` is a `Dict`, this returns a `Pair`.
        # We get the value entry using the `second` field of the `Pair`, which returns a `StreamHistory`,
        # of which we get the `complist` field.
        currentlist = first(A.list).second.complist
        X.complist != currentlist && error("all streams in StreamHistoryList must reference the same ComponentList")

        A.list[idx] = X
    end
end


function Base.getindex(A::StreamHistoryList, idx::String)
    return A.list[idx]
end


function Base.getindex(A::StreamHistoryList, idxs::Vector{String})
    res = StreamHistory[]
    for idx in idxs
        push!(res, A.list[idx])
    end
    return res
end


function Base.length(A::StreamHistoryList)
    return length(A.list)
end


function  Base.copy(A::Stream)
    return Stream(A.name, A.complist, A.comps, A.massflows, A.moleflows, A.totalmassflow, A.atomflows)
end


# Pretty printing for stream objects
function Base.show(io::IO, s::Stream)
    println(io, "Stream: ", s.name, " [Total mass flow: ", prettyround(s.totalmassflow), "]\n")
    println(io, "Component\tMass Flow\tMolar Flow")
    println(io, "-"^42)

    for (i, compname) in enumerate(s.comps)
        comp = s.complist[compname]
        println(io, " ", rpad(comp.name, 9), "\t", rpad(prettyround(s.massflows[i]), 9), "\t", prettyround(s.moleflows[i]))
    end

    println(io)
    println(io, "Atom\t\tFlow")
    println(io, "-"^20)
    atoms = collect(keys(s.atomflows)) 
    flows = collect(values(s.atomflows))

    for i in eachindex(atoms)
        println(io, " ", rpad(atoms[i],2), lpad(prettyround(flows[i]), 17))
    end
end


# Pretty printing for stream history objects
function Base.show(io::IO, s::StreamHistory)
    println(io, "StreamHistory: $(s.name)\n")
    println(io, "Components: \n")
    println(io, "-"^11)

    for compname in s.comps
        comp = s.complist[compname]
        println(io, " ", rpad(comp.name, 9))
    end

    println(io)
    println(io, "Data length: $(length(s.totalmassflow))")
    println(io, "Data starts:\t$(s.timestamps[begin])")
    println(io, "Data ends:\t$(s.timestamps[end])")
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


    @stream "mass" begin
        "Ethylene" --> 2.8053
        "Ethane" --> 27.06192
        "Hydrogen" --> 2.21738
    end syscomps "Test" sysstreams

    @stream "mole" begin 
        "Ethylene" --> 0.1
        "Ethane" --> 0.9
        "Hydrogen" --> 1.1
    end syscomps "Product" sysstreams 

"""
macro stream(flowtype::String, ex::Expr, complist::Symbol, name::String, streamlist::Symbol)      
    local comps = String[]
    local flows = Float64[]

    if lowercase(flowtype) == "mass"
        local ismoleflow = false
    elseif lowercase(flowtype) == "mole"
        local ismoleflow = true
    else
        error("flow basis specification must be mass or mole")
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
    return :($(esc(streamlist))[$name] = Stream($name, $(esc(complist)), $comps, $flows, $ismoleflow))
end


#----------------------------------------------------------------------------
# 
#----Utilities---------------------------------------------------------------
# 
#----------------------------------------------------------------------------


"""
    function copystream(list::StreamList, from::String, to::String)

Copy a stream in the stream list
"""
function copystream!(list::StreamList, from::String, to::String; factor=1.0)
    str = factor*list[from]
    list[to] = Stream(to, str.complist, str.comps, str.massflows, str.moleflows, str.totalmassflow, str.atomflows)

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
    str = list[from]
    list[to] = Stream(to, str.complist, str.comps, str.massflows, str.moleflows, str.totalmassflow, str.atomflows)
    delete!(list.list, from)

    return nothing
end


function renamestream(str::Stream, newname::String)
    return Stream(newname, str.complist, str.comps, str.massflows, str.moleflows, str.totalmassflow, str.atomflows)  
end


"""
    function copystreamhistory(list::StreamHistoryList, from::String, to::String)

Copy a StreamHistory in the StreamHistoryList
"""
function copystreamhistory!(list::StreamHistory, from::String, to::String; factor=1.0)
    str = factor*list[from]
    list[to] = StreamHistory(to, str.numdata, str.complist, str.comps, str.timestamps, str.massflowshistory, 
                             str.moleflowshistory, str.totalmassflowhistory, str.atomflowshistory)

    return nothing
end


"""
    function deletestreamhistory!(list::StreamHistoryList, from::String)

Delete a StreamHistory from the StreamHistoryList
"""
function deletestreamhistory!(list::StreamHistoryList, from::String)
    delete!(list.list, from)

    return nothing
end


"""
    function renamestream!(list::StreamHistoryList, from::String, to::String)

Rename a StreamHistory in the StreamHistoryList
"""
function renamestreamhistory!(list::StreamHistoryList, from::String, to::String)
    str = list[from]
    list[to] = StreamHistory(to, str.numdata, str.complist, str.comps, str.timestamps, str.massflowshistory, 
                             str.moleflowshistory, str.totalmassflowhistory, str.atomflowshistory)
    delete!(list.list, from)

    return nothing
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

    return StreamHistory(streamname, complist, comps, timestamps, transpose(flows), ismoleflow)
end








