#----------------------------------------------------------------------------
#
#----Definitions-------------------------------------------------------------
#
#----------------------------------------------------------------------------

struct BalanceBoundary
    name::String

    # Units included in the boundary
    unitops::UnitOpList
    included_units::Vector{String}

    # Number of historic data points
    numdata::Integer

    # Streams crossing the boundary limits
    # Calculate from inlets and outlets of included unit ops
    inlets::Vector{String}
    outlets::Vector{String}

    # Combined inlet and outlet streams
    total_in::Stream
    total_out::Stream

    # Total mass balance closure: out/in
    closure::TimeArray

    # Elemental closures: out/in
    atomclosures::TimeArray

    # Internal constructor to ensurte that all inlets and outlets have the same number of historic data points
    # and identical timestamps. It is only required here in case the user doesn't use the outer constructor
    function BalanceBoundary(name, units, included_units, numdata, inlets, outlets, total_in, total_out, closure, atomclosures)

        total_in.numdata != total_out.numdata && throw(DimensionMismatch("all in/outlets must must have similar history lengths."))

        !all(timestamp(total_in.massflows) .== timestamp(total_out.massflows)) && throw(DimensionMismatch("all in/outlets must must have identical timestamps."))

        new(name, units, included_units, numdata, inlets, outlets, total_in, total_out, closure, atomclosures)
    end
end

struct BoundaryList
    list::OrderedDict{String, BalanceBoundary}
end



#----------------------------------------------------------------------------
#
#----Constructors------------------------------------------------------------
#
#----------------------------------------------------------------------------


"""
    BalanceBoundary(name, unitops, included_units)

Constructor for a BalanceBoundary. Inputs are a UnitOpList and and array of names of `UnitOp`s in the list
that are inside the boundary. Inlet and outlet streams crossing the boundary are automatically calculated.

Will error if there are either no inlets or outlets.

Since not all atoms referenced in the streams will be present, closures for atoms not present will be indicated by
setting the values to -1.0
"""
function BalanceBoundary(name, unitops, included_units)
    # Get the streams that cross the boundary
    inlets, outlets, _ = boundarystreams(unitops, included_units)

    # Convert to list of names
    inletnames = String[]
    outletnames = String[]

    for inlet in inlets
        push!(inletnames, inlet.name)
    end
    for outlet in outlets
        push!(outletnames, outlet.name)
    end

    # Now add the inlets and outlets together to get one total inlet and one total outlet
    total_in = sum(inlets)
    total_out = sum(outlets)

    # Ensure that all inlets and outlets have the same number of historic data points at the same timestamps
    total_in.numdata != total_out.numdata && throw(DimensionMismatch("all in/outlets must must have similar history lengths."))
    !all(timestamp(total_in.massflows) .== timestamp(total_out.massflows)) && throw(DimensionMismatch("all in/outlets must must have identical timestamps."))

    # Calculate the closures
    closure = total_out.totalmassflow ./ total_in.totalmassflow
    rename!(closure, Symbol("Total Mass Closure"))
    numdata = length(closure)

    atomclosures = Vector{Dict{String, Float64}}(undef, numdata)

    for datum in 1:numdata
        _atomclosures = Dict{String, Float64}()
        for atom in keys(values(total_in.atomflows)[datum])
            in = values(total_in.atomflows)[datum][atom]
            out = values(total_out.atomflows)[datum][atom]
            _atomclosures[atom] = (out == 0.0) ? 0.0 : ((in == 0.0) ? -1.0 : out/in)
        end
        atomclosures[datum] = _atomclosures
    end
    atomclosures_ta = TimeArray(timestamp(closure), atomclosures, [Symbol("Elemental Closures")])
    BalanceBoundary(name, unitops, included_units, numdata, inletnames, outletnames, total_in, total_out, closure, atomclosures_ta)
end


"""
    Boundarylist()

Constructor for an empty boundary list. Boundaries are added when created via the @boundary macro
"""
function BoundaryList()
    l = OrderedDict{String, Stream}()
    return BoundaryList(l)
end


#----------------------------------------------------------------------------
#
#----Base overloads----------------------------------------------------------
#
#----------------------------------------------------------------------------


function Base.setindex!(A::BoundaryList, X::BalanceBoundary, idx::String)
    if length(A.list) == 0
        A.list[idx] = X
    else
        # Verify that all the streams reference the same UnitOpList.
        # Get the first item in the `list` field. As `list` is a `Dict`, this returns a `Pair`.
        # We get the value entry using the `second` field of the `Pair`, which returns a `Stream`,
        # of which we get the `complist` field.
        currentlist = first(A.list).second.unitops

        X.unitops != currentlist && throw(ArgumentError("all boundaries in BoundaryList must reference the same UnitOpList"))

        A.list[idx] = X
    end

    return nothing
end


function Base.getindex(A::BoundaryList, idx::String)
    return A.list[idx]
end


function Base.getindex(A::BoundaryList, idxs::Vector{String})
    res = BalanceBoundary[]
    for idx in idxs
        push!(res, A.list[idx])
    end
    return res
end


function Base.length(A::BoundaryList)
    return length(A.list)
end


# Pretty printing for BalanceBoundary objects
function Base.show(io::IO, b::BalanceBoundary)
    println(io, "Balance Boundary:\n")
    println(io, "Enclosed units: ", [u for u in b.included_units])
    println(io)
    println(io, "Closure:")
    pretty_table(io, b.closure, display_size=(14, -1))
    println(io)
    println(io, "Combined Feed Mass Flows:")
    pretty_table(io, b.total_in.massflows, display_size=(14, -1))
    println(io, "Combined Product Mass Flows:")
    pretty_table(io, b.total_out.massflows, display_size=(14, -1))
    println(io, "\nElemental closures (-1.0 if not present):")
    pretty_table(io, b.atomclosures, display_size=(14, -1))   
end


function  Base.show(io::IO, bl::BoundaryList)
    println(io, "Boundary list:")
    println(io)
    println(io, "Boundaries:")
    if length(bl.list) > 0
        for (name, _) in bl
            println(io, "  ", name)
        end 
    else
        println(io, "\tEmpty list")
    end
end


function Base.iterate(A::BoundaryList)
    return iterate(A.list)
end


function Base.iterate(A::BoundaryList, state)
    return iterate(A.list, state)
end

#----------------------------------------------------------------------------
#
#----Macros------------------------------------------------------------------
#
#----------------------------------------------------------------------------

"""

    @boundary begin
        unitops --> ["Reactor", "Membrane"]
    end "B1" fs

Create a boundary, that includes UnitOps "Reactor" amd "Membrane" from the UnitOpsList sysunitops and add it to BoundaryList sysboundaries.
Store it in fs.boundaries["B1"]

    @boundaryhist begin
        unitops --> ["RX101"]
    end "B2" fs

Create a boundary, that includes UnitOps "RX101" from the UnitOpsList sysunitops and add it to BoundaryList sysboundaries.
Store it in fs.boundaries["B2"]

"""
macro boundary(ex::Expr, name::String, fs::Symbol)      
    local units = String[]
    
    for line in ex.args
        match_comp = @capture(line, unitops --> [uops__])
        if match_comp
            for uop in uops
                push!(units, uop)
            end
        end
    end

    return :($(esc(fs)).boundaries[$name] = BalanceBoundary($name, $(esc(fs)).unitops, $units))
end


#----------------------------------------------------------------------------
# 
#----Utilities---------------------------------------------------------------
# 
#----------------------------------------------------------------------------


"""
    showdata(b::BalanceBoundary)

Returns a data table for a BalanceBoundary, as a String.

Example:
    print(showdata(b))
"""
function showdata(b::BalanceBoundary)
    invals = pretty_table(String, b.total_in.massflows)
    outvals = pretty_table(String, b.total_out.massflows)

    str = "Mass Balance Boundary:\n" * "-"^22 * "\n\nTotal Mass Flows In:\n" * invals * "\n\n" * "Total Mass Flows Out:\n" * outvals

    return str
end


"""
    inlets, outlets, internals = boundarystreams(units, included_units)

Find the streams that enter and leave the boundary that includes the specified unitops.
"""
function boundarystreams(unitops, included_units)
    #Figure out which streams cross the boundary:
    #    1. Take all the input stream from the units
    #    2. Subtract the ones that are product stream, as they will be internal streams
    #       Do the same for outlet streams, just the other way around
    #    3. Calculate the elemental closures

    inlets = Stream[]
    outlets = Stream[]
    internals = Stream[]

    # 1. Take all the input streams from the units 
    for unitname in included_units
        unit = unitops[unitname]

        for feedname in unit.inlets
            feed = unit.streams[feedname]
            push!(inlets, feed)
        end
        for prodname in unit.outlets
            prod = unit.streams[prodname]
            push!(outlets, prod)
        end
    end

    # 2a. Check which ones should be kept.
    # Don't delete now, or the check on outlets will miss the ones already deleted from inlets!
    keep_in = trues(length(inlets))
    keep_out = trues(length(outlets))
    for (i, stream) in enumerate(inlets)
        if stream in outlets
            keep_in[i] = false
            push!(internals, stream)
        end
    end
    for (i, stream) in enumerate(outlets)
        if stream in inlets
            keep_out[i] = false
        end
    end

    # 2b. Now delete the ones we don't want
    inlets = inlets[keep_in]
    outlets = outlets[keep_out]

    @argcheck length(inlets) > 0 "zero inlet streams to the boundary"
    @argcheck length(outlets) > 0 "zero outlet streams from the boundary"

    return inlets, outlets, internals
end

"""

    function deleteboundary!(fs, name)

Delete the specified bounary from the Flowsheet's BoundaryList.

"""
function deleteboundary!(fs, name)
    @argcheck fs isa Flowsheet "fs must be a Flowsheet"

    if name in keys(fs.boundaries.list)
        delete!(fs.boundaries.list, name);
    end
end


"""

    function deleteboundaries!(fs)

Delete all boundaries from the Flowsheet's BoundaryList.

"""
function deleteboundaries!(fs)
    @argcheck fs isa Flowsheet "fs must be a Flowsheet"

    for name in keys(fs.boundaries.list)
        delete!(fs.boundaries.list, name);
    end
end