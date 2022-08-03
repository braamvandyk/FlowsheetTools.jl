#----UnitOp------------------

struct UnitOp
    name::String

    # Connected streams
    inlets::Array{Stream, 1}
    outlets::Array{Stream, 1}
end


struct UnitOpHistory
    name::String

    # Number of historic data points
    numdata::Integer

    # Connected streams with history
    inlets::Array{StreamHistory, 1}
    outlets::Array{StreamHistory, 1}

    # Inner constructor to validate input data
    function UnitOpHistory(name, numdata, inlets, outlets)
        numdata = inlets[1].numdata
        timestamps = inlets[1].timestamps
        for inlet in inlets
            inlet.numdata != numdata && error("all streams must must have similar history lengths.")
            for j in eachindex(timestamps)
                timestamps[j] != inlet.timestamps[j] && error("all stream histories must have matching timestamps.")
            end
        end
        for outlet in outlets
            outlet.numdata != numdata && error("all streams must must have similar history lengths.")
            for j in eachindex(timestamps)
                timestamps[j] != inlet.timestamps[j] && error("all stream histories must have matching timestamps.")
            end
        end
        new(name, numdata, inlets, outlets)
    end
end


"""
    UnitOpHistory(name, inlets, outlets)

Constructor for a UnitOpHistory that defines the stream name and connected StreamHistory objects.
The internal field *numdata* is calculated from the connected streams.
"""
function UnitOpHistory(name, inlets, outlets)
    numdata = length(inlets[1].totalmassflow)

    UnitOpHistory(name, numdata, inlets, outlets)
end


# Pretty printing for UnitOp objects
function Base.show(io::IO, u::UnitOp)
    println(io, "Unit Operation: $(u.name)\n")
    println(io, "Feed streams:    ", [s.name for s in u.inlets])
    println(io)
    println(io, "Product streams: ", [s.name for s in u.outlets])
end


# Pretty printing for UnitOp objects
function Base.show(io::IO, u::UnitOpHistory)
    println(io, "Unit Operation: $(u.name)\n")
    println(io, "Feed streams:    ", [s.name for s in u.inlets])
    println(io)
    println(io, "Product streams: ", [s.name for s in u.outlets])
    println(io)
    println(io, "Data length:      $(length(u.numdata))")
end