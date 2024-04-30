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
function calccorrections_anchor(boundary::BalanceBoundary, anchor::String, customerror=nothing; totalweight=1.0, elementweight = 1.0, 
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
        !isnothing(customerror) && (totalerr += customerror(factors))
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
function calccorrections(boundary::BalanceBoundary, customerror=nothing; totalweight=1.0, elementweight = 1.0, 
    setelements = false, elementweights::Dict{String, Float64} = Dict{String, Float64}(), λ = 0.1)
    # Pull the streamlist of the first unit op in the list. Since this is from a UnitOpList,
    # all of the unit ops must have the same stream list
    streamlist = first(boundary.unitlist.list).second.streamlist  
     
    corrections = Dict{String, Float64}()

    inlets = boundary.inlets
    outlets = boundary.outlets
    streamnames = vcat(inlets, outlets)

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
        if !isnothing(customerror)
            calldict = Dict(zip(streamnames, factors))
            totalerr += customerror(calldict)
        end
        return totalerr
    end

    res = optimize(f, factors, LBFGS())
    # res = optimize(f, factors)
    optfactors = res.minimizer


    corrections = Dict(zip(vcat(inlets, outlets), optfactors))

    return corrections
end

# #=

function calccorrections(boundarylist::BoundaryList, customerror=nothing; totalweight=1.0, elementweight = 1.0, 
    setelements = false, elementweights::Dict{String, Float64} = Dict{String, Float64}(), λ = 0.1)
    # Pull the streamlist of the first unit op in first boundary in the list. Since this is from a UnitOpList,
    # all of the unit ops must have the same stream list
    firstboundary = first(values(boundarylist.list))
    streamlist = first(firstboundary.unitlist.list).second.streamlist  
    allinlets = Vector{Vector{String}}(undef, length(boundarylist))
    alloutlets = Vector{Vector{String}}(undef, length(boundarylist))
    allstreams  = String[]

    allins = Int64[]
    allouts = Int64[]

    for (i, boundary) in enumerate(boundarylist)
        allstreams = vcat(allstreams, boundary.second.inlets)
        allstreams = vcat(allstreams, boundary.second.outlets)
        allinlets[i] = boundary.second.inlets
        alloutlets[i] = boundary.second.outlets
        push!(allins, length(boundary.second.inlets))
        push!(allouts, length(boundary.second.outlets))
    end
    sort!(allstreams)
    unique!(allstreams)

    allfactors = ones(length(allstreams))

    numdata = streamlist[allinlets[1][1]].numdata


    function calcerr(boundarynum, allfactors)

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

        ins =  allins[boundarynum]
        outs = allouts[boundarynum]
        inlets = allinlets[boundarynum]
        outlets = alloutlets[boundarynum]
        
        factoridx_in = findall(x -> x in inlets, allstreams)
        factoridx_out = findall(x -> x in outlets, allstreams)

        factors = allfactors[vcat(factoridx_in, factoridx_out)]

        return f(factors)
    end

    function g(allfactors)
        totalerr = 0.0

        for boundarynum in 1:length(boundarylist)
            totalerr += calcerr(boundarynum, allfactors)
        end

        if !isnothing(customerror)
            calldict = Dict(zip(allstreams, allfactors))
            totalerr += customerror(calldict)
        end

        return totalerr
    end


    res = optimize(g, allfactors, LBFGS())
    optfactors = res.minimizer

    corrections = Dict(zip(allstreams, optfactors))

    return corrections
end

# =#


"""
    function closemb!(boundaries::BoundaryList; corrections::Dict{String, Float64})

Apply the mass balance reconciliation. The corrections can first calculated using `calccorrections`
and can then be applied to multiple boundaries using this function. If no corrections are passed,
`calccorrections` is automatically called.

Since balance boundaries are immutable, a new boundary instance is returned.
"""
function closemb!(boundaries::BoundaryList, corrections::Dict{String, Float64})
    # Pull the streamlist of the first unit op in the list. Since this is from a UnitOpList,
    # all of the unit ops must have the same stream list
    firstboundary = first(values(boundaries.list))
    unitlist = firstboundary.unitlist
    streamlist = first(unitlist.list).second.streamlist

    # Apply the stream corrections
    for stream in keys(corrections)
        streamlist[stream] *= corrections[stream]
    end

    #Recreate the boundaries
    for b in boundaries
        name = b.second.name
        units = b.second.units
        boundaries[name] = BalanceBoundary(name, unitlist, units)
    end

    return nothing
end


