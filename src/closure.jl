#----------------------------------------------------------------------------
#
#----Closure-----------------------------------------------------------------
#
#----------------------------------------------------------------------------


"""
    function calccorrections(boundary::BalanceBoundary; totalweight=1.0, elementweight=1.0)
    function calccorrections(boundary::BalanceBoundaryHistory; totalweight=1.0, elementweight=1.0)

Basic mass balance reconciliation. Error in total mass closure and average element closures are
weighted by *totalweight* and *elementweight* respectively and the squared weighted error is minimized.

Results are returned as a dict of streams and corrections to their flows.
"""
function calccorrections(boundary::BalanceBoundary; totalweight=1.0, elementweight=1.0)         
    # Pull the streamlist of the first unit op in the list. Since this is from a UnitOpList,
    # all of the unit ops must have the same stream list
    streamlist = first(boundary.unitlist.list).second.streamlist
    
    corrections = Dict{String, Float64}()
    
    inlets = boundary.inlets
    outlets = boundary.outlets

    ins = length(inlets)
    outs = length(outlets)
    
    factors = ones(ins + outs) # Correction factors for flows

    # Check for viable closure, i.e. there is at least some mass entering and some mass leaving.
    allzero = true
    for stream in inlets
        streamlist[stream].totalmassflow > 0.0 && (allzero = false)
    end
    allzero && throw(DomainError(x, "at least one inlet must have non-zero flow"))
   
    allzero = true
    for stream in inlets
        streamlist[stream].totalmassflow > 0.0 && (allzero = false)
    end
    allzero && throw(DomainError(x, "at least one inlet must have non-zero flow"))

    function f(factors)
        total_in = sum(factors[1:ins] .* streamlist[inlets])
        total_out = sum(factors[ins+1:end] .* streamlist[outlets])
        masserror = abs2(total_out.totalmassflow/total_in.totalmassflow - 1.0)

        atomclosures = Dict{String, Float64}()
        for atom in keys(total_in.atomflows)
            in = total_in.atomflows[atom]
            out = total_out.atomflows[atom]
            atomclosures[atom] = out/in
        end
        numatoms = length(atomclosures)
        atomerror = abs2(sum(values(atomclosures)) - numatoms)
        return totalweight*masserror + elementweight*atomerror
    end


    res = optimize(f, factors)
    optfactors = res.minimizer

    i = 1
    for stream in inlets
        corrections[stream] = optfactors[i]
        i += 1
    end
    for stream in outlets
        corrections[stream] = optfactors[i]
        i += 1
    end

    return corrections
end


"""

    function closemb(boundary::BalanceBoundary, [corrections::Dict{String, Float64}])
    function closemb(boundary::BalanceBoundaryHistory, [corrections = Dict{String, Float64}])

Apply the mass balance reconciliation. The corrections can first calculated using `calccorrections`
and can then be applied to multiple boundaries using this function. If no corrections are passed,
`calccorrections` is automatically called.

Since balance boundaries are immutable, a new boundary instance is returned.
"""
function closemb(boundary::BalanceBoundary, corrections=nothing; totalweight=1.0, elementweight=1.0)
    isnothing(corrections) && (corrections = calccorrections(boundary; totalweight, elementweight))

    # Pull the streamlist of the first unit op in the list. Since this is from a UnitOpList,
    # all of the unit ops must have the same stream list
    streamlist = first(boundary.unitlist.list).second.streamlist  

    # Apply the stream corrections
    for stream in keys(corrections)
        streamlist[stream] *= corrections[stream]
    end

    #Recreate the boundary
    newboundary = BalanceBoundary(boundary.unitlist, boundary.units)

    return newboundary
end



function calccorrections(boundary::BalanceBoundaryHistory; totalweight=1.0, elementweight=1.0)         
    # Pull the streamlist of the first unit op in the list. Since this is from a UnitOpList,
    # all of the unit ops must have the same stream list
    streamlist = first(boundary.unitlist.list).second.streamlist  
     
    corrections = Dict{String, Float64}()

    inlets = boundary.inlets
    outlets = boundary.outlets

    ins = length(inlets)
    outs = length(outlets)

    factors = ones(ins + outs)

    numdata = streamlist[inlets[1]].numdata
   
    function f(factors)
        # Since inlets and outlets are arrays of StreamHistory, summing them produces StreamHistory objects
        total_in = sum(factors[1:ins] .* streamlist[inlets])
        total_out = sum(factors[ins+1:end] .* streamlist[outlets])
        
        masserror = 0.0
        atomerror = 0.0
        for datum = 1:numdata
            masserror += abs2(total_out.totalmassflow[datum]/total_in.totalmassflow[datum] - 1.0)
            
            atomclosures = Dict{String, Float64}()
            for atom in keys(total_in.atomflows[datum])
                in = total_in.atomflows[datum][atom]
                out = total_out.atomflows[datum][atom]
                atomclosures[atom] = out/in
            end
            numatoms = length(atomclosures)
            atomerror += abs2(sum(values(atomclosures)) - numatoms)
        end
        return totalweight*masserror + elementweight*atomerror
    end

    res = optimize(f, factors)
    optfactors = res.minimizer

    i = 1
    for stream in inlets
        corrections[stream] = optfactors[i]
        i += 1
    end
    for stream in outlets
        corrections[stream] = optfactors[i]
        i += 1
    end

    return corrections
end


function closemb(boundary::BalanceBoundaryHistory, corrections=nothing; totalweight=1.0, elementweight=1.0)
    isnothing(corrections) && (corrections = calccorrections(boundary; totalweight, elementweight))

    # Pull the streamlist of the first unit op in the list. Since this is from a UnitOpList,
    # all of the unit ops must have the same stream list
    streamlist = first(boundary.unitlist.list).second.streamlist  

    # Apply the stream corrections
    for stream in keys(corrections)
        streamlist[stream] *= corrections[stream]
    end

    #Recreate the boundary
    newboundary = BalanceBoundaryHistory(boundary.unitlist, boundary.units)

    return newboundary
end