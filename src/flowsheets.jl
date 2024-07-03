#----------------------------------------------------------------------------
#
#----Definitions-------------------------------------------------------------
#
#----------------------------------------------------------------------------

struct Flowsheet
    comps::ComponentList
    streams::StreamList
    unitops::UnitOpList
    boundaries::BoundaryList

    rununits::Vector{String}
    runorder::Vector{Int}
end


function (fs::Flowsheet)(neworder = nothing; showoutput=true)
    if isnothing(neworder)
        for i in fs.execution_order
            showoutput && println(i, "  ", fs.unitops[i])
            fs.unitlist[fs.unitops[i]]()
        end
    else
        for i in neworder
            u = fs.unitops[i]
            showoutput && println(i, "  ", fs.unitops[i])
            u()
        end
    end

    return nothing
end


#----------------------------------------------------------------------------
#
#----Constructors------------------------------------------------------------
#
#----------------------------------------------------------------------------


function Flowsheet()
    complist = ComponentList()
    streamlist = StreamList()
    unitlist = UnitOpList()
    boundarylist = BoundaryList()

    return Flowsheet(complist, streamlist, unitlist, boundarylist, String[], Int[])
end


#----------------------------------------------------------------------------
# 
#----Base overloads----------------------------------------------------------
# 
#----------------------------------------------------------------------------


function Base.show(io::IO, fs::Flowsheet)
    d = OrderedDict{Int, String}()
    for i in 1:length(fs.runorder)
        d[fs.execution_order[i]] = fs.unitops[i]
    end
    sort!(d)

    println(io, "Flowsheet:")
    println(io)
    println(io, "Execution Order:")
    for (k, v) in d
        println(k, "  ", v)
    end
end


#----------------------------------------------------------------------------
# 
#----Utilities---------------------------------------------------------------
# 
#----------------------------------------------------------------------------


"""
    addtorun!(fs::Flowsheet, u::String)


Add a single unit operation, or an array of unit operations to the execution list of a Flowsheet object.
Default execution order is set as the order in which the unit operations are added.
"""
function addtorun!(fs::Flowsheet, u::String)
    # Only add if it is not already there.
    if isnothing(findfirst(x -> x == u, fs.rununits))
        if haskey(fs.unitops.list, u)
            push!(fs.rununits, u)
            push!(fs.runorder, length(fs.runorder) + 1)
        else
            throw(ArgumentError("unit operation $u is not defined in list of unit operations for $fs"))
        end
    end

    return nothing
end

"""
    removeunitop!(fs::Flowsheet, u::String)
    removeunitop!(fs::Flowsheet, ::Vector{String})

Removed the specified unit operation(s) from the Flowsheet object, as well as from its execution order.
"""
function removeunitop!(fs::Flowsheet, u::String)
    index = findfirst(x -> x == u, fs.unitops)
    #  Ignore if not there
    if !isnothing(index)
        popat!(fs.unitops, index)
        for i in fs.execution_order
            #  Remove from the execution order
            if i == index
                popat!(fs.execution_order, i)
            end
            #  All entries > the one to remove, must move down one position
            if i > index
                fs.execution_order[i] -= 1
            end
        end
    end
end


function removeunitop!(fs::Flowsheet, ::Vector{String})
    for u in us
        index = findfirst(x -> x == u, fs.unitops)
        #  Ignore if not there
        if !isnothing(index)
            popat!(fs.unitops, index)
            for i in fs.execution_order
                #  Remove from the execution order
                if i == index
                    popat!(fs.execution_order, i)
                end
                #  All entries > the one to remove, must move down one position
                if i > index
                    fs.execution_order[i] -= 1
                end
            end
        end
    end
end


"""
    setorder!(fs::Flowsheet, neworder)

Change the execution order of the unit operations in a Flowsheet object.
The lenght of the execution order array does not have to match the number of units. Indicidual units can be executed multiple times, or not at all.
"""
function setorder!(fs::Flowsheet, neworder)
    numold = length(fs.execution_order)
    numnew = length(neworder)

    if numold == numnew
        for i in 1:numold
            fs.execution_order[i] = neworder[i]
        end
    elseif numold < numnew
        for i in 1:numold
            fs.execution_order[i] = neworder[i]
        end            
        for j = 1:(numnew - numold)
            push!(fs.execution_order, neworder[numold+j])
        end
    else
        for i in 1:numold
            fs.execution_order[i] = neworder[i]
        end            
        for j = 1:(numold - numnew)
            pop!(fs.execution_order)
        end
    end

    return nothing
end


"""
    generateBFD(fs::Flowsheet, filename = nothing; displaybfd=true)

Generate a block flow diagram of the flowsheet (using Mermaid.js).
If a filename is specified, the diagram will be downloaded to that file, else a filename generated by Mermaid.js will be used.
The file type can be either `.svg` or `.png`. Specifying other extensions to the filename will result in an exception.
If no filename is specified, the file will be an SVG image.

if `displaybfd` is true, the diagram will also be displayed in the default manner for the environment, using `display()`.
"""
function generateBFD(fs::Flowsheet, filename = nothing; displaybfd=true)
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

    inlets, outlets, internals = boundarystreams(fs.unitops, fs.rununits)
    
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
        for unitname in fs.rununits
            unit = fs.unitops[unitname]
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
        for unitname in fs.rununits
            unit = fs.unitops[unitname]
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
        for unitname in fs.rununits
            unit = fs.unitops[unitname]
            if outlet in unit.outlets
                push!(streams, (replace(unitname, r"\s+" =>s""), replace(outlet, r"\s+" =>s""), " "))
                break
            end
        end
    end

    mermaidlines = "%%{init: {'securityLevel': 'loose', 'theme': 'base', 'themeVariables': { 'fontSize': '14', 'darkMode': 'true', 'primaryColor': '#1e90ff', 'lineColor': '#1e90ff'}}}%%\ngraph LR\n"
    # mermaidlines = "{'fontSize': '14', 'darkMode': 'true', 'primaryColor': '#1e90ff', 'lineColor': '#1e90ff'}}}%%\ngraph LR\n"    
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

    # To generate the Mermaid diagram, encode the instructions in base 64,
    # generate a url and download the file from mermaid.ink.

    # First check the file extension, if filename specified
    if !isnothing(filename)
        name, ext = splitext(filename)
        ext = lowercase(ext[2:end])
        if ext ∉ ["svg", "png"]
            throw(ArgumentError("Can only generate .svg or .png files."))
        end
        (ext == "png") && (ext = "img")
    else
        ext = "svg"
    end
    
    io = IOBuffer()
    iob64_encode = Base64EncodePipe(io)
    write(iob64_encode, mermaidlines)
    close(iob64_encode)
    mmstr = String(take!(io))
    mmurl = "http://mermaid.ink/" * ext * "/" * mmstr
    
    if !isnothing(filename)
        Downloads.download(mmurl, filename)
    else
        filename = Downloads.download(mmurl)
    end

    if ext == "svg"
        f = open(filename, "r")
        img = read(f, String)
        close(f)
        displaybfd && display("image/svg+xml", img)
    else
        f = open(filename, "r")
        img = read(f)
        close(f)
        displaybfd && display("image/png", img)
    end
  
    # Return the filename
    return filename
end



"""

    componentnames(fs)

Return the list of names of components in the included `ComponentList` in the `Flowsheet`, fs.

# Examples
```julia
julia> componentnames(fs)
```

"""
function componentnames(fs)
    return collect(keys(fs.comps.list))
end