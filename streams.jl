#----Stream------------------

struct Stream
    name::String

    comps::Array{Component, 1}
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
    
    comps::Array{Component, 1}

    timestamps::Array{DateTime, 1}
    # Specifiy mass flows, calculate mole flows
    massflows::Array{Float64, 2}
    moleflows::Array{Float64, 2}
    totalmassflow::Array{Float64, 1}

    # Molar flow of atoms, calculated from above
    atomflows::Array{Dict{String, Float64}, 1}
end


"""
    Stream(name, comps, massflows)

Constructor for a stream that defines the stream name and component flows.
The mass flows are specified and a molar composition and atomic molar flows calculated.
"""
function Stream(name, comps, massflows)
    numcomps = length(comps)
    moleflows = zeros(numcomps)
    atomflows = Dict{String, Float64}()
    totalmassflow = 0.0
    
    for (i, comp) in enumerate(comps)
        totalmassflow += massflows[i]
        moleflows[i] = massflows[i]/comp.Mr
        
        for (j, atom) in enumerate(comp.atoms)
            if atom in keys(atomflows)
                atomflows[atom] += moleflows[i]*comp.counts[j]
            else
                atomflows[atom] = moleflows[i]*comp.counts[j]
            end
        end
    end

    Stream(name, comps, massflows, moleflows, totalmassflow, atomflows)
end


"""
    StreamHistory(name, comps, timestamps, massflowshistory)

Constructor for a stream history object that defines the stream name and component flows
for various past measurements.
The mass flows are specified and a molar composition and atomic molar flows calculated.
The mass flow history is passed as a matrix where each column is a datum
"""
function StreamHistory(name, comps, timestamps, massflowshistory)
    numcomps = length(comps)
    numdata = size(massflowshistory, 2)

    length(timestamps) != numdata && error("length mismatch between data and timestamps.")

    moleflowshistory = zeros(numcomps, numdata)
    atomflowshistory = Array{Dict{String, Float64}, 1}(undef, numdata)
    totalmassflowhistory = zeros(numdata)
    
    for datum in 1:numdata
        atomflows = Dict{String, Float64}()

        for (i, comp) in enumerate(comps)
            totalmassflowhistory[datum] += massflowshistory[i, datum]
            moleflowshistory[i, datum] = massflowshistory[i, datum]/comp.Mr
            
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

    StreamHistory(name, numdata, comps, timestamps, massflowshistory, moleflowshistory, totalmassflowhistory, atomflowshistory)
end


"""
    Base.+(a::Stream, b::Stream)

Extend the addition operator to add to streams to each other - a mixer.
It is assumed that the streams will have different components in arbitrary order.
"""
function Base.:+(a::Stream, b::Stream)
    comps = copy(a.comps)
    massflows = copy(a.massflows)

    numcomps = length(comps)

    compmap = indexin(b.comps, comps)
    for (i, comp) in enumerate(b.comps)
        if isnothing(compmap[i])
            # Add the component
            push!(comps, b.comps[i])
            push!(massflows, b.massflows[i])
        else
            # Add the flow to the existing component
            massflows[compmap[i]] += b.massflows[i]
        end
    end

    Stream(a.name * "-" * b.name, comps, massflows)
end


"""
    Base.+(a::StreamHistory, b::StreamHistory)

Extend the addition operator to add to streams histories to each other - a mixer.
It is assumed that the streams histories will have different components in arbitrary order.
"""
function Base.:+(a::StreamHistory, b::StreamHistory)
    # Check that the data lengths are the same.
    # The constructor already checks that the timestamp length matches the data length.
    a.numdata != b.numdata && throw(DimensionMismatch("length of a not equal to size of b"))

    # Check that the timestamps are identical!
    for i in eachindex(a.timestamps)
        a.timestamps[i] != b.timestamps[i] && error("timestamp values do not match at entry $i")
    end


    comps = copy(a.comps)
    massflows = copy(a.massflows)

    numcomps = length(comps)
    numdata = size(massflows, 2)

    compmap = indexin(b.comps, comps)
    for (i, comp) in enumerate(b.comps)
        if isnothing(compmap[i])
            # Add the component
            push!(comps, b.comps[i])
            vcat(massflows, b.massflows[i, :]')
        else
            # Add the flow to the existing component
            massflows[compmap[i], :] .+= b.massflows[i, :]
        end
    end

    StreamHistory(a.name * "-" * b.name, comps, a.timestamps, massflows)
end


"""
    Base.*(a::T, b::Stream) where T <: Real

Extend the multiplication operator to scale a stream's flows by a scalar value.
Used in mass balance reconciliations to apply flow corrections.
"""
function Base.:*(a::T, b::Stream) where T <: Real
    Stream(b.name, b.comps, a .* b.massflows)
end


"""
    Base.*(a::T, b::StreamHistory) where T <: Real

Extend the multiplication operator to scale a stream history's flows by a scalar value.
Used in mass balance reconciliations to apply flow corrections.
"""
function Base.:*(a::T, b::StreamHistory) where T <: Real
    StreamHistory(b.name, b.comps, b.timestamps, a .* b.massflows)
end


"""
    function copystream(s::Stream, name::String)

Copy a stream. Can also be used to rename a stream.
"""
function copystream(s::Stream, name::String)
    str = Stream(name, s.comps, s.massflows, s.moleflows, s.totalmassflow, s.atomflows)
end


"""
    function copystream(s::Stream, name::String)

Copy a stream. Can also be used to rename a stream.
"""
function copystreamhistory(s::StreamHistory, name::String)
    str = StreamHistory(name, s.numdata, s.comps, s.timestamps, s.massflows, s.moleflows, s.totalmassflow, s.atomflows)
end


# Pretty printing for stream objects
function Base.show(io::IO, s::Stream)
    println(io, "Stream: ", s.name, " [Total mass flow: ", prettyround(s.totalmassflow), "]\n")
    println(io, "Component\tMass Flow\tMolar Flow")
    println(io, "-"^42)
    for (i, comp) in enumerate(s.comps)
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
    for comp in s.comps
        println(io, " ", rpad(comp.name, 9))
    end
    println(io)
    println(io, "Data length: $(length(s.totalmassflow))")
end


# Parsing an expression to extract components and massflows for use in macro @stream
function parse_stream(ex, name)
    comps = Component[]
    massflows = Float64[]
    
    for line in ex.args
        match_comp = @capture(line, comp_sym_ --> flow_)
        if match_comp
            comp = eval(comp_sym)
            massflow = eval(flow)
            if any(x -> x == comp, comps)
                i = findfirst(x -> x == comp, comps)
                massflows[i] += massflow
            else
                push!(comps, comp)
                push!(massflows, massflow)
            end
        end
    end
    return :(stream = Stream($name, $comps, $massflows))
end


"""
    @stream begin
        ethylene --> 2.0
        hydrogen --> 6.2
    end "Feed"

    @stream begin
        ethylene --> ethylene.Mr
        hydrogen --> 2.0*hydrogen.Mr
    end "Feed"

Defines a Stream() with the specified name and component mass flows.
The flows may be expressions.
"""
macro stream(ex::Expr, name::String)      
    return parse_stream(ex, name)
end


function readstreamhistory(streamname, comps)
    df = CSV.read(joinpath("streamhistories", streamname * ".csv"), DataFrame, dateformat="yyyy/mm/dd HH:MM")
    massflows = zeros(size(df, 1), length(comps))
    
    for (i, comp) in enumerate(comps)
        massflows[:, i] = df[:, comp.name]
    end

    return StreamHistory(streamname, comps, df.TimeStamp, transpose(massflows))
end
