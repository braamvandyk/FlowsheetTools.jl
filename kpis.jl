#----KPIs--------------------

# KPIs are now defined as functions that take a BalanceBoundary as input.
# All mass / mole / atom flows in and out are in the BalanceBoundary, 
# with mass and elemental closures already calculated.

"""
    conversion(b::BalanceBoundary, c::Component)

Calculate the fractional conversion of the component (c) over the balance boundary (b)
"""
function conversion(b::BalanceBoundary, c::Component)
    # Find the component in the feed
    i = findfirst(x->x==c, b.total_in.comps)
    if isnothing(i)
        feed = 0.0
    else
        feed = b.total_in.massflows[i]
    end

    # Find the component in the product
    i = findfirst(x->x==c, b.total_out.comps)
    if isnothing(i)
        prod = 0
    else
        prod = b.total_out.massflows[i]
    end

    return feed > 0.0 ? (feed - prod)/feed : 0.0
end


"""
    function selectivity(b, reactant, product)

Calculate the molar selectivity of the converted reactant to the product, over the balance boundary (b).
"""
function selectivity(b, reactant, product)
    # Find the reactant inflow
    i = findfirst(x->x==reactant, b.total_in.comps)
    if isnothing(i)
        r_in = 0.0
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