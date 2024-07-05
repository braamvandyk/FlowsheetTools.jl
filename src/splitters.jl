#----------------------------------------------------------------------------
#
#----Splitters--------------------------------------------------------------
#
#----------------------------------------------------------------------------


"""

    flowsplitter!(fs, outlets::Vector{String}, inlets::Vector{String}, params)

Calculation for flowsplitter UnitOps. Splits combined feed stream according to params - an iterable with 
split fractions to (n - 1) streams. The last stream gets the balance.

Example:

    @unitop begin
        inlets --> ["Product"]
        outlets --> ["Product1", "Product2"]
        calc --> flowsplitter!
        params --> [0.5]
    end "ProductSplitter" fs

"""
function flowsplitter!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)
    numouts = length(outlets)
    sumfracs = sum(params)

    @argcheck length(params) == (numouts - 1) "incorrect number of fractions - must be one less than outlets"
    @argcheck 0.0 <= sumfracs <= 1.0 "invalid split fractions specified"

    totalin = sum(streamlist[inlets])
    
    for (i, outname) in enumerate(outlets)
        if i < numouts
            frac = params[i]
        else
            frac = 1.0 - sumfracs
        end

        streamlist[outname] = renamestream(frac * totalin, outname)
    end
    
    return nothing
end


"""
    
    componentplitter!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)

Calculation for component splitter UnitOps. Splits combined feed stream according to params - a nested Dict with 
split fractions for each component to (n - 1) streams. The last stream gets the balance.

Example:

    @unitop begin
    inlets --> ["Product1"]
    outlets --> ["Product1a", "Product1b"]
    calc --> componentplitter!
    params --> Dict([
        "Hydrogen" => Dict(["Product1a" => 0.5]),
        "Ethane" => Dict(["Product1b" => 0.3])
    ])
    end "ComponentSplitter" fs

"""
function componentplitter!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)
    complist = first(streamlist).second.comps
    numcomps = length(complist)
    numouts = length(outlets)
    totalin = sum(streamlist[inlets])
    fractions = fill(-1.0, (numcomps, numouts))
    compnames = Array{String}(undef, numcomps)

    # Loop through the specified component splits and assign the specified fractions to their streams
    # Any fraction not specified will remain as -1.0, which will be processed in the next step
    for (compindex, comp) in enumerate(complist)
        compname = comp.first
        compnames[compindex] = compname
        if compname in keys(params)
            #  This component has specified splits
            splits = params[compname]
            @argcheck length(splits) == numouts - 1 "incorrect number of fractions specified for component $compname"
            for split in splits
                stream = split.first
                @argcheck stream in outlets "specified stream not in list of outlet streams"
                streamindex = findfirst(str -> str == stream, outlets)
                fraction = split.second
                @argcheck 0.0 <= fraction <= 1.0 "invalid fraction specified for $compname in stream @stream"
                fractions[compindex, streamindex] = fraction
            end
        end
    end
    
    #  Now we look for unspecified fractions and set those.
    for row in 1:size(fractions)[1]
        sumfracs = sum(fractions[row, :])
        if sumfracs == -numouts
            # No splits specified for this component -> everything goes to last stream
            fractions[row, :] .= 0.0
            fractions[row, end] = 1.0
        else
            # We already know from the previous step that n - 1 streams are assigned for
            # each component with specified splits, or the assert would have failed.
            # This means the only unspecified component still has a value of -1.0. This
            # fraction must be the 1 - sum(specified fractions)
            for j in 1:numouts
                if fractions[row, j] == -1.0
                    fractions[row, j] = 1.0 - (sumfracs + 1.0)
                    break
                end
            end
        end

    end

    total_massflows = values(totalin.massflows)
    timestamps = timestamp(totalin.massflows)
    for streamindex in 1:numouts
        massflows = copy(total_massflows)
        for datum in 1:size(total_massflows)[1]
            massflows[datum, :] .*= fractions[:, streamindex]
        end
        streamlist[outlets[streamindex]] = Stream(outlets[streamindex], complist, compnames, timestamps, massflows)
    end

    return nothing
end