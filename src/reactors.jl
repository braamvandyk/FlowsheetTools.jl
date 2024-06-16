# TODO Add a check of reaction stoichiometry to the constructor

#----------------------------------------------------------------------------
#
#----Definitions-------------------------------------------------------------
#
#----------------------------------------------------------------------------


struct Reaction
    complist:: ComponentList

    reactants::Array{String}
    products::Array{String}
    reactcoeffs::Array{Float64}
    prodcoeffs::Array{Float64}
    
    targetcomp::String
    targetconversion::Float64

    targetindex::Int
end


#----------------------------------------------------------------------------
#
#----Pretty Printing---------------------------------------------------------
#
#----------------------------------------------------------------------------


function Base.show(io::IO, r::Reaction)
    print(io, "Reaction:\t")
    for (i, reactant) in enumerate(r.reactants)
        print(io, r.reactcoeffs[i], " ", reactant)
        if i < length(r.reactants)
            print(io, " + ")
        end
    end
    print(io, " --> ")
    for (i, product) in enumerate(r.products)
        print(io, r.prodcoeffs[i], " ", product)
        if i < length(r.products)
            print(io, " + ")
        end
    end

    println(io)
    println(io, "Conversion of $(r.targetcomp) is $(r.targetconversion)")
end



#----------------------------------------------------------------------------
#
#----Constructors------------------------------------------------------------
#
#----------------------------------------------------------------------------


"""
    function Reaction(reactants, products, reactcoeffs, prodcoeffs, targetcomp, targetconversion)

Constructor for Reaction object. Calculates the index of the target reactant in the reactants vector and stores it for efficiency in downstream calculations.
"""
function Reaction(complist, reactants, products, reactcoeffs, prodcoeffs, targetcomp, targetconversion)
    targetindex = findfirst(isequal(targetcomp), reactants)

    # Do stoichiometry check
    atoms_in = Dict{String, Float64}()
    atoms_out = Dict{String, Float64}()

    for (i, reactant) in enumerate(reactants)
        for (j, atom) in enumerate(complist[reactant].atoms)
            atoms_in[atom] = get(atoms_in, atom, 0.0) + reactcoeffs[i]*complist[reactant].counts[j]
        end
    end
    for (i, product) in enumerate(products)
        for (j, atom) in enumerate(complist[product].atoms)
            atoms_out[atom] = get(atoms_out, atom, 0.0) + prodcoeffs[i]*complist[product].counts[j]
        end
    end

    !(atoms_in == atoms_out) && throw(ArgumentError("Reaction stoichiometry check failed."))
    
    return Reaction(complist, reactants, products, reactcoeffs, prodcoeffs, targetcomp, targetconversion, targetindex)
end


#----------------------------------------------------------------------------
#
#----Reactors----------------------------------------------------------------
#
#----------------------------------------------------------------------------


"""
    stoichiometric_reactor!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)

A simple reactor block with specified reaction stoichiometry and fractional conversions for limiting components.
All reactions are assumed to happen in parallel and therefore all fractional conversions are relative to the total feed.     

"""
function stoichiometric_reactor!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)
    @argcheck length(outlets) == 1 "only one outlet stream allowed."
    numreactions = length(params)

    targetcomps = Array{String}(undef, numreactions)
    
    for (i, reaction) in enumerate(params)
        @argcheck  isa(reaction, Reaction) "parameters passed to the reactor must be of type Reaction"
        targetcomps[i] = reaction.targetcomp
    end
    unique!(targetcomps)
    
    sumconversions = zeros(length(targetcomps))
    for reaction in params
        index = findfirst(==(reaction.targetcomp), targetcomps)
        sumconversions[index] += reaction.targetconversion
    end     
 
    @argcheck all(sumconversions .<= 1.0) "total conversion for $(targetcomps[sumconversions .> 1.0]) is greater than 1.0"
    
    # Combine all feed streams
    feed = sum(streamlist[inlets])
    moleflows = feed.moleflows

    # Calculate all the fractional conversions for all components (not just limiting ones)

    # For each reaction, get the delta moleflow for each species
    # Add these deltas up for all of the reactions and then apply to the feed
    # Save the result in the product

    # Change in molar flow of each possible species
    # Remember that moleflows is a TimeArray
    delta_moleflows = Dict{String, Vector{Float64}}()
    for reaction in params
        targetcomp = reaction.targetcomp
        targetconversion = reaction.targetconversion
        targetindex = reaction.targetindex # index of the target reactant in the list of reactants for this reaction
        targetcoeff = reaction.reactcoeffs[targetindex]

        for (i, reactant) in enumerate(reaction.reactants)
            # Remember that moleflows is a TimeArray, so access the values via values()
            if reactant in keys(delta_moleflows)
                delta_moleflows[reactant] .-= values(moleflows[targetcomp]) .* targetconversion .* reaction.reactcoeffs[i] ./ targetcoeff
            else
                delta_moleflows[reactant] = -1.0 .* values(moleflows[targetcomp]) .* targetconversion .* reaction.reactcoeffs[i] ./ targetcoeff
            end
        end
        for (i, product) in enumerate(reaction.products)
            if product in keys(delta_moleflows)
                delta_moleflows[product] .+= values(moleflows[targetcomp]) .* targetconversion .* reaction.prodcoeffs[i] ./ targetcoeff
            else
                delta_moleflows[product] = values(moleflows[targetcomp]) .* targetconversion .* reaction.prodcoeffs[i] ./ targetcoeff
            end
        end
    end
        
    # Assign to the exit stream
    timestamps = timestamp(moleflows)
    comps = string.(colnames(moleflows))
    vals = zeros(length(timestamps), length(comps))
    # Calculate the exit stream composition
    # Add the values in delta_moleflows to the values in moleflows and create an output stream
    for component in keys(delta_moleflows)
        idx = findfirst(isequal(component), comps)
        vals[:,idx] = values(moleflows[component]) + delta_moleflows[component]
    end
    
    # Replace outlet stream in streamlist. Set ismoleflow to true.
    streamlist[outlets[1]] = Stream(outlets[1], feed.complist, comps, timestamps, vals, true)  
    
    return nothing
end