#----------------------------------------------------------------------------
#
#----Definitions-------------------------------------------------------------
#
#----------------------------------------------------------------------------

struct Flowsheet
    unitlist::UnitOpList
    units::Vector{String}

    order::Vector{Int}
end

function (fs::Flowsheet)(neworder = nothing)
    if isnothing(neworder)
        for i in fs.order
            fs.unitlist[fs.units[i]]()
        end
    else
        for i in neworder
            u = fs.units[i]
            u()
        end
    end

    return nothing
end


#----------------------------------------------------------------------------
# 
#----Utilities---------------------------------------------------------------
# 
#----------------------------------------------------------------------------


function addunitop!(fs::Flowsheet, u::String)
    if haskey(fs.unitlist.list, u)
        push!(fs.units, u)
        push!(fs.order, length(fs.order) + 1)
    else
        error("unitop $u not defined in list $(fs.unitlist)")
    end

    return nothing
end


function addunitop!(fs::Flowsheet, us::Vector{String})
    for u in us
        if haskey(fs.unitlist.list, u)
            push!(fs.units, u)
            push!(fs.order, length(fs.order) + 1)
        else
            error("unitop $u not defined in list $(fs.unitlist)")
        end
    end

    return nothing
end


function setorder!(fs::Flowsheet, neworder)
    numold = length(fs.order)
    numnew = length(neworder)

    if numold == numnew
        for i in 1:numold
            fs.order[i] = neworder[i]
        end
    elseif numold < numnew
        for i in 1:numold
            fs.order[i] = neworder[i]
        end            
        for j = 1:(numnew - numold)
            push!(fs.order, neworder[numold+j])
        end
    else
        for i in 1:numold
            fs.order[i] = neworder[i]
        end            
        for j = 1:(numold - numnew)
            pop!(fs.order)
        end
    end

    return nothing
end


function generateBFD(fs::Flowsheet)
#=
    Assign each unit op an id = its name (without whitespace)    
    
    1. Call boundarystreams() to return feed, internal and product streams
    2. or each stream, return its source and destination
        for feed/product streams, the source/destination has the same id, but append "port", or similar

    
    

    Iterate over all inlet streams and add <inlet id>( ) -->|inlet stream name| <block id>[block name]
    Iterate over all internal streams and add <source id> -->|internal stream name| <id>[block name]   (the source blocks should already have names)
    Iterate over all outlet streams andd add <id> -->|outlet stream name| <outlet id>( )


    graph LR
    feed1( ) -->|stream 1| block1[B1]
    block1 -->|stream 2| block2[B2]
    block1 -->|stream 3| block3[B3]
        
=#

    inlets, outlets, internals = boundarystreams(fs.unitlist, fs.units)
    
    # Convert to list of names
    inletnames = String[]
    outletnames = String[]
    internalnames = String[]
    
    for inlet in inlets
        push!(inletnames, inlet.name)
    end
    for outlet in outlets
        push!(outletnames, outlet.name)
    end
    for internal in internals
        push!(internalnames, internal.name)
    end


    streams = Tuple{String, String, String}[]
    for inlet in inletnames
        # Find the block that this stream connects to
        for unitname in fs.units
            unit = fs.unitlist[unitname]
            if inlet in unit.inlets
                push!(streams, (" ", replace(inlet, r"\s+" =>s""), replace(unitname, r"\s+" =>s"")))
                break
            end
        end
    end

    for internal in internalnames
        foundsource = false
        foundsink = false
        thissource = ""
        thissink = ""

        # Find the blocks that this stream connects to
        for unitname in fs.units
            unit = fs.unitlist[unitname]
            if internal in unit.outlets
                thissource = unitname
                foundsource = true
            end
            if internal in unit.inlets
                thissink = unitname
                foundsink = true
            end
            (foundsink & foundsource) && break           
        end
        push!(streams, (replace(thissource, r"\s+" =>s""), replace(internal, r"\s+" =>s""), replace(thissink, r"\s+" =>s"")))
    end

    for outlet in outletnames
        # Find the block that this stream connects to
        for unitname in fs.units
            unit = fs.unitlist[unitname]
            if outlet in unit.outlets
                push!(streams, (replace(unitname, r"\s+" =>s""), replace(outlet, r"\s+" =>s""), " "))
                break
            end
        end
    end

    mermaidlines = "%%{init: {'securityLevel': 'loose', 'theme': 'default', 'themeVariables': { 'fontSize': '14'}}}%%\ngraph LR\n"
    for stream in streams
        if stream[1] == " " # Feed stream
            mermaidline = "$(stream[2])([ ]) -->|$(stream[2])| $(stream[3])[$(stream[3])]"
        elseif stream[3] == " " # Product stream
            mermaidline = "$(stream[1])[$(stream[1])] -->|$(stream[2])| $(stream[2])([ ])"
        else
            mermaidline = "$(stream[1])[$(stream[1])] -->|$(stream[2])| $(stream[3])[$(stream[3])]"
        end

        mermaidlines *= "$mermaidline\n"
    end

    return Diagram(:mermaid, mermaidlines)

end