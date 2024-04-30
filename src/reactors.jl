struct Reaction
    complist::ComponentList
    reactants::Array{String}
    products::Array{String}
    reactcoeffs::Array{Float64}
    prodcoeffs::Array{Float64}
    
    targetcomp::String
    targetconversion::Float64
end




"""
    stoichiometric_reactor!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)

A simple reactor block with specified reaction stoichiometry and fractional conversions for limiting components.
All reactions are assumed to happen in parallel and therefore all fractional conversions are relative to the total feed.     

"""
function stoichiometric_reactor!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)
    # @argcheck length(outlets) == 1 "only one outlet stream allowed."
    # numreactions = length(params)

    # targetcomps = Array{String}(length(params))
    # for (i, reaction) in enumerate(params)
    #     @argcheck  isa(reaction, Reaction) "parameters passed to the reactor must be of type Reaction"
    #     targetcomps[i] = reaction.targetcomp
    # end
    # unique!(targetcomps)
    # @argcheck lenght(targetcomps) == length(params) "target components not unique for all reactions"

    

    # # Combine all feed streams
    # feed = sum(streamlist[inlets])
    # moleflows = feed.moleflows

    # # Calculate the exit stream composition
    # for reaction in params
    #     reacted = moleflows[reaction.targetcomp] .* reaction.targetconversion
    #     moleflows[reaction.targetcomp] .= moleflows[reaction.targetcomp] - reacted

    #     # for each reaction, find the reactants and their coefficients
    #     # for each non-target reactant, ratio their coefficient with the target comp and get the moles recreated
    #     # Check if the resulting moles left is not negative - error is so?
    #     # Now calculate the product moles and add to what is in the outlet stream
    # end




    # # Assign to the exit stream
    # name = streamlist[outlets[1]].name     
    # tempstream = sum(streamlist[inlets])


    
    
    return nothing
end