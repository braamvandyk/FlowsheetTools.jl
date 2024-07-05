#----------------------------------------------------------------------------
#
#----Closure-----------------------------------------------------------------
#
#----------------------------------------------------------------------------


"""

    calccorrections(fs; customerror=nothing, anchor = nothing, totalweight=1.0, elementweight = 1.0, 
        setelements = false, elementweights::Dict{String, Float64} = Dict{String, Float64}(), λ = 0.1)

Basic mass balance reconciliation. Error in total mass closure and overall element closures are
weighted by *totalweight* and *elementweight* respectively and the squared weighted error is minimized.
If *setelements* is true, the dictionary of weights for each element is applied, rather than value of *elementweights*.

One stream can be specified as an anchor stream, that will not be modified during reconciliation - correction factor of 1.0.
This stream is identified by supplying it's name in the kwarg `anchor`.

Regularisation, as per ridge regression is done with a weight, λ = 0.1 by default.

Results are returned as a dict of stream names and the correction factors for their flows.

"""
function calccorrections(fs; customerror=nothing, anchor = nothing, totalweight=1.0, elementweight = 1.0, 
    setelements = false, elementweights::Dict{String, Float64} = Dict{String, Float64}(), λ = 0.1)


    # Places to store names on inlets and outlets for each boundary
    allinlets = Vector{Vector{String}}(undef, length(fs.boundaries))
    alloutlets = Vector{Vector{String}}(undef, length(fs.boundaries))
    
    # Places to store counts of inlets and outlets for each boundary
    allins = Int64[]
    allouts = Int64[]
    
    # Place to store the names of all non-anchor streams
    # This will be used to index into the correction factors vector
    allstreams  = String[]

    for (i, boundary) in enumerate(fs.boundaries)
        # boundary is a name => value pair from the boundary list, so use boundary.second to get the actual boundary object
        allstreams = vcat(allstreams, boundary.second.inlets)
        allstreams = vcat(allstreams, boundary.second.outlets)
        allinlets[i] = boundary.second.inlets
        alloutlets[i] = boundary.second.outlets
        push!(allins, length(boundary.second.inlets))
        push!(allouts, length(boundary.second.outlets))
    end
    sort!(allstreams)
    unique!(allstreams)

    # Now exclude the anchor, if any
    if !isnothing(anchor)
        anchorindex = searchsorted(allstreams, anchor)
        deleteat!(allstreams, anchorindex)
    end
    # Start off with all factors = 1, so no corrections
    allfactors = ones(length(allstreams))

    numdata = fs.streams[allinlets[1][1]].numdata


    function boundaryerror(boundarynum, allfactors)
    # Extract the streams involved with this boundary, look up the relevant factors
    # and then calculate the error from this boundary.

        function thiserror(factors)
            # This function receives the factors for all streams, with the anchor set to one,
            # if present. First inlets, then outlets.

            # Since inlets and outlets are arrays of Stream, summing them produces Stream objects
            total_in = sum(factors[1:ins] .* fs.streams[inlets])
            total_out = sum(factors[ins+1:end] .* fs.streams[outlets])
            
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

        # Get this boundary's inlets and outlets
        # ins and outs are the counts of inlets and outlets
        ins =  allins[boundarynum]
        outs = allouts[boundarynum]
        # inlets and outlets are the names of the inlet and outlet streams
        inlets = allinlets[boundarynum]
        outlets = alloutlets[boundarynum]

        # Now get the factors for each stream [inlets outlets] and remember to set the anchor
        # to one, if it is in this boundary
        factors = Array{Float64}(undef, ins + outs)
        for (index, inlet) in enumerate(inlets)
            if inlet == anchor
                factors[index] = 1.0
            else
                allstreams_index = findfirst(isequal(inlet), allstreams)
                factors[index] = allfactors[allstreams_index]
            end
        end
        for (index, outlet) in enumerate(outlets)
            if outlet == anchor
                factors[ins + index] = 1.0
            else
                allstreams_index = findfirst(isequal(outlet), allstreams)
                factors[ins + index] = allfactors[allstreams_index]
            end
        end


        return thiserror(factors)
    end

    function loss(allfactors)
        totalerr = 0.0

        for boundarynum in 1:length(fs.boundaries)
            totalerr += boundaryerror(boundarynum, allfactors)
        end

        if !isnothing(customerror)
            calldict = Dict(zip(allstreams, allfactors))
            totalerr += customerror(calldict)
        end

        return totalerr
    end


    res = optimize(loss, allfactors, BFGS())
    optfactors = res.minimizer

    corrections = Dict(zip(allstreams, optfactors))

    return corrections
end


"""

    closemb!(fs; corrections::Dict{String, Float64})

Apply the mass balance reconciliation correction factors. The corrections can first calculated using
`calccorrections` and can then be applied to multiple boundaries using this function. 

"""
function closemb!(fs, corrections::Dict{String, Float64})
    # Apply the stream corrections
    for stream in keys(corrections)
        fs.streams[stream] *= corrections[stream]
    end

    #Recreate the boundaries
    for b in fs.boundaries
        name = b.second.name
        units = b.second.included_units
        fs.boundaries[name] = BalanceBoundary(name, fs.unitops, units)
    end

    return nothing
end


