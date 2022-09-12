#----------------------------------------------------------------------------
#
#----Definitions-------------------------------------------------------------
#
#----------------------------------------------------------------------------

struct BalanceBoundary
    # Units included in the boundary
    unitlist::UnitOpList
    units::Vector{String}

    # Streams crossing the boundary limits
    # Calculate from inlets and outlets of included unit ops
    inlets::Vector{String}
    outlets::Vector{String}

    # Combined inlet and outlet streams
    total_in::Stream
    total_out::Stream

    # Total mass balance closure: (out - in)/in
    closure::Float64

    # Elemental closures, Dict{String, Float64}
    atomclosures::Dict{String, Float64}
end


struct BalanceBoundaryHistory
    # Units included in the boundary
    unitlist::UnitOpHistoryList
    units::Vector{String}

    # Number of historic data points
    numdata::Integer

    # Streams crossing the boundary limits
    # Calculate from inlets and outlets of included unit ops
    inlets::Vector{String}
    outlets::Vector{String}

    # Combined inlet and outlet streams
    total_in::StreamHistory
    total_out::StreamHistory

    # Total mass balance closure: (out - in)/in
    closure::Vector{Float64}

    # Elemental closures, Dict{String, Float64}
    atomclosures::Vector{Dict{String, Float64}}


    function BalanceBoundaryHistory(unitlist, units, numdata, inlets, outlets, total_in, total_out, closure, atomclosures)
        total_in.numdata != total_out.numdata && error("All in/outlets must must have similar history lengths.")
        new(unitlist, units, numdata, inlets, outlets, total_in, total_out, closure, atomclosures)
    end
end


#----------------------------------------------------------------------------
#
#----Constructors------------------------------------------------------------
#
#----------------------------------------------------------------------------


"""
    BalanceBoundary(unitlist, units)

Constructor for a BalanceBoundary. Inputs are a UnitOpList and and array of names of UnitOps in the list
that are inside the boundary. Inlet and outlet streams crossing the boundary are automatically calculated.

Fill error if there are either no inlets or outlets.
"""
function BalanceBoundary(unitlist::UnitOpList, units::Vector{String})
    # Get the streams that cross the boundary
    inlets, outlets = boundarystreams(unitlist, units)

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

    # Calculate the closures
    closure = total_out.totalmassflow/total_in.totalmassflow

    atomclosures = Dict{String, Float64}()
    for atom in keys(total_in.atomflows)
        in = total_in.atomflows[atom]
        out = total_out.atomflows[atom]
        atomclosures[atom] = out/in
    end

    BalanceBoundary(unitlist, units, inletnames, outletnames, total_in, total_out, closure, atomclosures)
end


"""
    BalanceBoundaryHistory(unitlist, units)

Constructor for a BalanceBoundaryHistory. Inputs are a UnitOpHistoryList and and array of names of `UnitOpHistory`s in the list
that are inside the boundary. Inlet and outlet streams crossing the boundary are automatically calculated.

Fill error if there are either no inlets or outlets.
"""
function BalanceBoundaryHistory(unitlist::UnitOpHistoryList, units::Vector{String})
    # Get the streams that cross the boundary
    inlets, outlets = boundarystreams(unitlist, units)

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

    # Calculate the closures
    closure = total_out.totalmassflow./total_in.totalmassflow
    numdata = length(closure)

    atomclosurehistory = Vector{Dict{String, Float64}}(undef, numdata)


    for datum in 1:numdata
        atomclosures = Dict{String, Float64}()
        for atom in keys(total_in.atomflows[datum])
            in = total_in.atomflows[datum][atom]
            out = total_out.atomflows[datum][atom]
            atomclosures[atom] = out/in
        end
        atomclosurehistory[datum] = atomclosures
    end

    BalanceBoundaryHistory(unitlist, units, numdata, inletnames, outletnames, total_in, total_out, closure, atomclosurehistory)
end


#----------------------------------------------------------------------------
#
#----Base overloads----------------------------------------------------------
#
#----------------------------------------------------------------------------


# Pretty printing for BalanceBoundary objects
function Base.show(io::IO, b::BalanceBoundary)
    println(io, "Balance Boundary:\n")
    println(io, "Enclosed units: ", [u for u in b.units])
    println(io, "Closure: ", prettyround(b.closure))
    println(io)
    print(io, "Combined Feed ")
    println(io, b.total_in)
    print(io, "Combined Product ")
    println(io, b.total_out)
    println(io, "\nElemental closures:")
    println(io, "="^19)

    atoms = collect(keys(b.atomclosures)) 
    closures = collect(values(b.atomclosures))
    for i in eachindex(atoms)
        println(io, " ", rpad(atoms[i],2), lpad(prettyround(closures[i]), 14))
    end    
end


# Pretty printing for BalanceBoundary objects
function Base.show(io::IO, b::BalanceBoundaryHistory)
    println(io, "Balance Boundary:\n")
    println(io, "Enclosed units: ", [u for u in b.units])
    println(io)
    println(io, "Data length:\t", b.numdata)
    println(io, "Data starts:\t", b.total_in.timestamps[begin])
    println(io, "Data ends:\t", b.total_in.timestamps[begin])
end


#----------------------------------------------------------------------------
#
#----Macros------------------------------------------------------------------
#
#----------------------------------------------------------------------------

"""
    @boundary begin
        unitops --> ["Reactor", "Membrane"]
    end b sysunitops

Create a boundary, b, that includes UnitOps "Reactor" amd "Membrane" from the UnitOpsList sysunitops.
"""
macro boundary(ex::Expr, name::Symbol, unitoplist::Symbol)      
    local unitops = String[]
    
    for line in ex.args
        match_comp = @capture(line, unitops --> [ous__])
        if match_comp
            for ou in ous
                push!(unitops, ou)
            end
        end
    end

    return :($(esc(name)) = BalanceBoundary($(esc(unitoplist)), $unitops))
end


"""
    @boundaryhist begin
        unitops --> ["RX101"]
    end bh histops

Create a boundary, bh, that includes UnitOps "RX101" from the UnitOpsList histops.
"""
macro boundaryhist(ex::Expr, name::Symbol, unitoplist::Symbol)      
    local unitops = String[]
    
    for line in ex.args
        match_comp = @capture(line, unitops --> [ous__])
        if match_comp
            for ou in ous
                push!(unitops, ou)
            end
        end
    end

    return :($(esc(name)) = BalanceBoundaryHistory($(esc(unitoplist)), $unitops))
end


#----------------------------------------------------------------------------
# 
#----Utilities---------------------------------------------------------------
# 
#----------------------------------------------------------------------------


"""
    showdata(bh::BalanceBoundaryHistory)

Returns a data table for a BalanceBoundaryHistory, as a String.

Example:
    print(showdata(bh))
"""
function showdata(bh::BalanceBoundaryHistory)
    titles_in = copy(bh.total_in.comps)
    titles_in .= "In:" .* titles_in
    titles_out = copy(vcat(bh.total_out.comps))
    titles_out .= "Out:" .* titles_out
    titles = vcat("Timestamp", titles_in, titles_out)

    data_in = hcat(bh.total_in.timestamps, bh.total_in.massflows')
    data_out = bh.total_out.massflows'
    data = hcat(data_in, data_out)

    str = pretty_table(String, data, header=titles)

    return str
end


"""
    inlets, outlets = boundarystreams(unitlist, units)

Find the streams that enter and leave the boundary that includes the specified unitops.
"""
function boundarystreams(unitlist::UnitOpList, units::Vector{String})
    #Figure out which streams cross the boundary:
    #    1. Take all the input stream from the units
    #    2. Subtract the ones that are product stream, as they will be internal streams
    #       Do the same for outlet streams, just the other way around
    #    3. Calculate the elemental closures

    inlets = Stream[]
    outlets = Stream[]
    internals = Stream[]

    # 1. Take all the input stream from the units 
    for unitname in units
        unit = unitlist[unitname]

        for feedname in unit.inlets
            feed = unit.streamlist[feedname]
            push!(inlets, feed)
        end
        for prodname in unit.outlets
            prod = unit.streamlist[prodname]
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
            push!(internals, stream) # only need this in this loop, as any internal stream is some block's inlet
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

    @assert length(inlets) > 0 "zero inlet streams to the boundary"
    @assert length(outlets) > 0 "zero outlet streams from the boundary"

    return inlets, outlets, internals
end


function boundarystreams(unitlist::UnitOpHistoryList, units::Vector{String})
    #Figure out which streams cross the boundary:
    #    1. Take all the input stream from the units
    #    2. Subtract the ones that are product stream, as they will be internal streams
    #       Do the same for outlet streams, just the other way around
    #    3. Calculate the elemental closures

    inlets = StreamHistory[]
    outlets = StreamHistory[]

    # 1. Take all the input stream from the units 
    for unitname in units
        unit = unitlist[unitname]

        for feedname in unit.inlets
            feed = unit.streamlist[feedname]
            push!(inlets, feed)
        end
        for prodname in unit.outlets
            prod = unit.streamlist[prodname]
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

    @assert length(inlets) > 0 "zero inlet streams to the boundary"
    @assert length(outlets) > 0 "zero outlet streams from the boundary"

    return inlets, outlets
end