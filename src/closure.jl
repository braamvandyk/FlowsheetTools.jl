#----------------------------------------------------------------------------
#
#----Closure-----------------------------------------------------------------
#
#----------------------------------------------------------------------------



"""
    function calccorrections_anchor(boundary::BalanceBoundary, anchor::String; totalweight=1.0, elementweight = 1.0, 
        setelements = false, elementweights::Dict{String, Float64} = Dict{String, Float64}) 

Basic mass balance reconciliation. Error in total mass closure and average element closures are
weighted by *totalweight* and *elementweight* respectively and the squared weighted error is minimized.

Results are returned as a dict of streams and corrections to their flows.
"""
function calccorrections_anchor(boundary::BalanceBoundary, anchor::String; totalweight=1.0, elementweight = 1.0, 
    setelements = false, elementweights::Dict{String, Float64} = Dict{String, Float64}()) 
    # Pull the streamlist of the first unit op in the list. Since this is from a UnitOpList,
    # all of the unit ops must have the same stream list
    streamlist = first(boundary.unitlist.list).second.streamlist 

    corrections = Dict{String, Float64}()

    inlets = boundary.inlets
    outlets = boundary.outlets

    ins = length(inlets)
    outs = length(outlets)

    if anchor in inlets
        anchorinlet = true
    elseif anchor in outlets
        anchorinlet = false
    else
        throw(ArgumentError("anchor $anchor not in inlets or outlets"))
    end

    if anchorinlet
        anchorindex = findfirst(isequal(anchor), inlets)
    else
        anchorindex = findfirst(isequal(anchor), outlets) + ins
    end
    

    factors = ones(ins + outs - 1) # No factor for the anchor stream

    numdata = streamlist[inlets[1]].numdata
   
    function f(_factors)
        factors = vcat(_factors[1:anchorindex-1], 1.0, _factors[anchorindex:end])

        # Since inlets and outlets are arrays of Stream, summing them produces Stream objects
        total_in = sum(factors[1:ins] .* streamlist[inlets])
        total_out = sum(factors[ins+1:end] .* streamlist[outlets])
        
        masserrors = (total_out.totalmassflow ./ total_in.totalmassflow) .- 1.0
        masserror = sum(abs2, values(masserrors))
        
        atomerror = 0.0
        for datum = 1:numdata           
            for atom in keys(values(total_in.atomflows)[datum])
                inflow = values(total_in.atomflows)[datum][atom]
                if inflow > 0.0
                    outflow = values(total_out.atomflows)[datum][atom]
                    if setelements
                        aweight = !(atom ∈ keys(elementweights)) ? 0.0 : elementweights[atom]
                    else
                        aweight = elementweight
                    end
                    atomerror +=  aweight * abs2(outflow/inflow - 1.0)
                end
            end
        end
        totalerr = totalweight*masserror + atomerror
        return totalerr
    end

    res = optimize(f, factors, BFGS())
    # res = optimize(f, factors)
    _optfactors = res.minimizer
    optfactors = vcat(_optfactors[1:anchorindex-1], 1.0, _optfactors[anchorindex:end])

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
    function calccorrections(boundary::BalanceBoundary; totalweight=1.0, elementweight = 1.0, 
        setelements = false, elementweights::Dict{String, Float64}= Dict{String, Float64}, λ = 0.1)

Basic mass balance reconciliation. Error in total mass closure and average element closures are
weighted by *totalweight* and *elementweight* respectively and the squared weighted error is minimized.

Results are returned as a dict of streams and corrections to their flows.
"""
function calccorrections(boundary::BalanceBoundary; totalweight=1.0, elementweight = 1.0, 
    setelements = false, elementweights::Dict{String, Float64} = Dict{String, Float64}(), λ = 0.1)
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
        # Since inlets and outlets are arrays of Stream, summing them produces Stream objects
        total_in = sum(factors[1:ins] .* streamlist[inlets])
        total_out = sum(factors[ins+1:end] .* streamlist[outlets])
        
        masserrors = (total_out.totalmassflow ./ total_in.totalmassflow) .- 1.0
        masserror = sum(abs2, values(masserrors))
        
        atomerror = 0.0
        for datum = 1:numdata           
            for atom in keys(values(total_in.atomflows)[datum])
                inflow = values(total_in.atomflows)[datum][atom]
                if inflow > 0.0
                    if setelements
                        aweight = !(atom ∈ keys(elementweights)) ? 0.0 : elementweights[atom]
                    else
                        aweight = aweight = elementweight
                    end
                    outflow = values(total_out.atomflows)[datum][atom]
                    atomerror += aweight * abs2(outflow/inflow - 1.0)
                end
            end
        end
        totalerr = totalweight*masserror + atomerror + λ*sum(abs2, 1.0 .- factors)
        return totalerr
    end

    res = optimize(f, factors, LBFGS())
    # res = optimize(f, factors)
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
    function closemb_anchor(boundary::BalanceBoundary; corrections=nothing, anchor=nothing, totalweight=1.0, elementweight=1.0)

Apply the mass balance reconciliation. The corrections can first calculated using `calccorrections`
and can then be applied to multiple boundaries using this function. If no corrections are passed,
`calccorrections` is automatically called.

Since balance boundaries are immutable, a new boundary instance is returned.
"""
function closemb(boundary::BalanceBoundary, corrections::Dict{String, Float64})
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
