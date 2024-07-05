#----KPIs--------------------

# KPIs are defined as functions that take a BalanceBoundary as input.
# All mass / mole / atom flows in and out are in the BalanceBoundary, 
# with mass and elemental closures already calculated.


"""

    conversion(b::BalanceBoundary, component::String)

Calculate the fractional conversion of the component over the balance boundary.
Results are returned for each timestamp in the streams.

Example

    conversion(fs.boundaries["B1"], "Ethane")

"""
function conversion(b::BalanceBoundary, component::String)
    # Find the component in the feed
    if component in keys(b.total_in.comps.list)
        feeds = values(b.total_in.massflows[Symbol(component)])
    else
        feeds = zeros(length(b.total_in.massflows))
    end

    # Find the component in the product
    if component in keys(b.total_in.comps.list)
        prods = values(b.total_out.massflows[Symbol(component)])
    else
        prods = zeros(length(b.total_out.massflows))
    end

    conversions = zeros(length(b.total_in.massflows))
    for (i, (feed, prod)) in enumerate(zip(feeds, prods))
        if feed > 0.0
            conversions[i] = (feed - prod)/feed
        else
            conversions[i] = 0.0
        end
    end

    if length(conversions) == 1
        return conversions[1]
    else
        timestamps = timestamp(b.total_in.massflows)
        conversions_ta = TimeArray(timestamps, conversions, [Symbol("Conversion: " * component)])
        return conversions_ta
    end
end


"""

    molar_selectivity(b, reactant, product)

Calculate the molar selectivity of the converted reactant to the product, over the balance boundary (b).
Since the calculation does not know about actual reactions, it calculates Δproduct / Δreactant.
Results are returned for each timestamp in the streams.

Example
    molar_selectivity(fs.boundaries["B1"], "Ethylene", "Ethane")

"""
function molar_selectivity(b::BalanceBoundary, reactant::String, product::String)
    # Find the reactant inflow
    if reactant in keys(b.total_in.comps.list)
        r_ins = values(b.total_in.moleflows[Symbol(reactant)])
    else
        r_ins = zeros(length(b.total_in.moleflows))
    end

    # Find the reactant outflow
    if reactant in keys(b.total_out.comps.list)
        r_outs = values(b.total_out.moleflows[Symbol(reactant)])
    else
        r_ins = zeros(length(b.total_out.moleflows))
    end

    # Find the product inflow
    if product in keys(b.total_in.comps.list)
        p_ins = values(b.total_in.moleflows[Symbol(product)])
    else
        p_ins = zeros(length(b.total_in.moleflows))
    end

    # Find the product outflow
    if product in keys(b.total_out.comps.list)
        p_outs = values(b.total_out.moleflows[Symbol(product)])
    else
        p_ins = zeros(length(b.total_out.moleflows))
    end

    selectivities = zeros(length(b.total_in.massflows))
    for (i, (r_in, r_out, p_in, p_out)) in enumerate(zip(r_ins, r_outs, p_ins, p_outs))
        if r_in > 0.0
            selectivities[i] = (p_out - p_in)/(r_in - r_out)
        else
            selectivities[i] = 0.0
        end
    end

    if length(selectivities) == 1
        return selectivities[1]
    else
        timestamps = timestamp(b.total_in.massflows)
        selectivities_ta = TimeArray(timestamps, selectivities, [Symbol("Molar selectivity: " * reactant * " -> " * product)])
        return selectivities_ta
    end

    return selectivities
end
