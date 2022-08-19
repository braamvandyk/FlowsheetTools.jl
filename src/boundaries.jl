#----------------------------------------------------------------------------
#
#----Definitions-------------------------------------------------------------
#
#----------------------------------------------------------------------------

struct BalanceBoundary
    # Units included in the boundary
    unitlist::UnitOpList
    units::Array{String, 1}

    # Streams crossing the boundary limits
    # Calculate from inlets and outlets of included unit ops
    inlets::Array{String, 1}
    outlets::Array{String, 1}

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
    units::Array{String, 1}

    # Number of historic data points
    numdata::Integer

    # Streams crossing the boundary limits
    # Calculate from inlets and outlets of included unit ops
    inlets::Array{String, 1}
    outlets::Array{String, 1}

    # Combined inlet and outlet streams
    total_in::StreamHistory
    total_out::StreamHistory

    # Total mass balance closure: (out - in)/in
    closure::Array{Float64, 1}

    # Elemental closures, Dict{String, Float64}
    atomclosures::Array{Dict{String, Float64}, 1}


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


function BalanceBoundary(unitlist, units)

#Figure out which streams cross the boundary:
#    1. Take all the input stream from the units
#    2. Subtract the ones that are product stream, as they will be internal streams
#       Do the same for outlet streams, just the other way around
#    3. Calculate the elemental closures

    inlets = Stream[]
    outlets = Stream[]

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

    # 2c. And convert to list of names
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


function BalanceBoundaryHistory(unitlist::UnitOpHistoryList, units::Vector{String})
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

    # 2c. And convert to list of names
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

    atomclosurehistory = Array{Dict{String, Float64}, 1}(undef, numdata)


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
    println(io, "Data length:    ", b.numdata)
end