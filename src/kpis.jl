#----KPIs--------------------

# KPIs are now defined as functions that take a BalanceBoundary as input.
# All mass / mole / atom flows in and out are in the BalanceBoundary, 
# with mass and elemental closures already calculated.


"""
    conversion(b::BalanceBoundary, component::String)
    conversion(b::BalanceBoundaryHistory, component::String)

Calculate the fractional conversion of the component over the balance boundary
"""
function conversion(b::BalanceBoundary, component::String)
    # Find the component in the feed
    i = findfirst(x->x==component, b.total_in.comps)
    if isnothing(i)
        feed = 0.0
    else
        feed = b.total_in.massflows[i]
    end

    # Find the component in the product
    i = findfirst(x->x==component, b.total_out.comps)
    if isnothing(i)
        prod = 0
    else
        prod = b.total_out.massflows[i]
    end

    return feed > 0.0 ? (feed - prod)/feed : 0.0
end


function conversion(b::BalanceBoundaryHistory, component::String) 
    # Find the component in the feed
    i = findfirst(x->x==component, b.total_in.comps)
    
    # Conversion zero if nothing in the feed
    if isnothing(i) 
        results = zeros(b.total_in.numdata)
        return results
    end

    # Find the component in the product
    j = findfirst(x->x==component, b.total_out.comps)

    # If nothing left in the product, then return ones
    if isnothing(j)
        results = ones(b.total_in.numdata)
        return results
    end
    
    # Create place to put the results.
    results = zeros(b.total_in.numdata)

    for datum = 1:b.total_in.numdata
        feed = b.total_in.massflows[i, datum]
        prod = b.total_out.massflows[j, datum]
        results[datum] = feed > 0.0 ? (feed - prod)/feed : 0.0
    end

    return results
end


"""
    function selectivity(b, reactant, product)
    function selectivity(b::BalanceBoundaryHistory, reactant::String, product::String)

Calculate the molar selectivity of the converted reactant to the product, over the balance boundary (b).
"""
function selectivity(b::BalanceBoundary, reactant::String, product::String)
    # Find the reactant inflow
    i = findfirst(x->x==reactant, b.total_in.comps)
    if isnothing(i)
        return 0.0
    else
        r_in = b.total_in.moleflows[i]
    end

    # Find the reactant outflow
    i = findfirst(x->x==reactant, b.total_out.comps)
    if isnothing(i)
        r_out = 0.0
    else
        r_out = b.total_out.moleflows[i]
    end

    # Find the product inflow
    i = findfirst(x->x==product, b.total_in.comps)
    if isnothing(i)
        p_in = 0
    else
        p_in = b.total_in.moleflows[i]
    end

    # Find the product outflow
    i = findfirst(x->x==product, b.total_out.comps)
    if isnothing(i)
        p_out = 0.0
    else
        p_out = b.total_out.moleflows[i]
    end

    return r_in > 0.0 ? (p_out - p_in)/(r_in - r_out) : 0.0
end


function selectivity(b::BalanceBoundaryHistory, reactant::String, product::String)
    results = zeros(b.total_in.numdata)

    # Find the reactant in the feed
    i = findfirst(x->x==reactant, b.total_in.comps)
    # Selectivity zero if nothing in the feed
    isnothing(i) && (return results)
    r_in = b.total_in.moleflows[i, :]
    
    # Find the reactant in the product
    i = findfirst(x->x==reactant, b.total_out.comps)
    if isnothing(i)
        r_out = zeros(b.total_out.numdata)
    else
        r_out = b.total_out.moleflows[i, :]
    end

    # Find the product in the feed
    i = findfirst(x->x==product, b.total_in.comps)
    if isnothing(i)
        p_in = zeros(b.total_in.numdata)
    else
        p_in = b.total_in.moleflows[i, :]
    end
    
    # Find the product in the product
    i = findfirst(x->x==product, b.total_out.comps)
    if isnothing(i)
        p_out = zeros(b.total_out.numdata)
    else
        p_out = b.total_out.moleflows[i, :]
    end

    for datum = 1:b.total_in.numdata
        results[datum] = r_in[datum] > 0.0 ? (p_out[datum] - p_in[datum])/(r_in[datum] - r_out[datum]) : 0.0
    end

    return results
end