using FlowsheetTools
using Graphs
using GraphMakie, CairoMakie



fs = Flowsheet()
_ = readcomponentlist!(fs, "components", ["Ethylene", "Ethane", "Hydrogen"])

# @stream mole begin
#     "Hydrogen" --> 1.1
# end "H2" fs

# @stream mole begin
#     "Ethylene" --> 0.1
#     "Ethane" --> 0.9
# end "C2" fs

# @stream mole begin
#     "Ethylene" --> 0.0
#     "Ethane" --> 1.0
#     "Hydrogen" --> 1.0
# end "Product" fs

# addemptystream!(fs, "Mixed")
# addemptystream!(fs, "RXOutlet")
# addemptystream!(fs, "Product")
# addemptystream!(fs, "Recycle")

# @unitop begin
#     inlets --> ["H2", "C2", "Recycle"]
#     outlets --> ["Mixed"]
#     calc --> mixer!
# end "Mixer" fs

# @unitop begin
#     inlets --> ["Mixed"]
#     outlets --> ["RXOutlet"]
# end "Reactor" fs

# @unitop begin
#     inlets --> ["RXOutlet"]
#     outlets --> ["Product", "Recycle"]
#     calc --> flowsplitter!
#     params --> [0.8]
# end "RecycleSplitter" fs

# fs()

# generateBFD(fs, "./teartest.svg")

@stream mole begin
    "Ethylene" --> 0.0
    "Ethane" --> 1.0
    "Hydrogen" --> 1.0
end "1" fs
addemptystream!(fs, "2")
addemptystream!(fs, "3")
addemptystream!(fs, "4")
addemptystream!(fs, "5")
addemptystream!(fs, "6")
addemptystream!(fs, "7")
addemptystream!(fs, "8")

@unitop begin
    inlets --> ["7"]
    outlets --> ["1"]
end "A" fs

@unitop begin
    inlets --> ["1", "8"]
    outlets --> ["2"]
end "B" fs

@unitop begin
    inlets --> ["2", "5"]
    outlets --> ["3", "6", "7"]
end "C" fs

@unitop begin
    inlets --> ["3"]
    outlets --> ["4", "8"]
end "D" fs

@unitop begin
    inlets --> ["4", "6"]
    outlets --> ["5"]
end "E" fs



# --------------------------------
# 
# This code will replace the same in flowsheets.jl once it works
# 
# --------------------------------


function incidence_matrix(fs::Flowsheet)
    # Create a directed graph with the unit operations as nodes
    numnodes = length(fs.unitops.list)
    nodenames = collect(keys(fs.unitops.list))
    numedges = length(fs.streams.list)
    edgenames = collect(keys(fs.streams.list))

    inc_matrix = zeros(Int, numnodes, numedges)

    for (j, unitpair) in enumerate(fs.unitops )
        unit = unitpair.second #the UnitOp object
        for i in 1:length(unit.inlets)
            streamname = unit.inlets[i]
            idx = findall(==(streamname), edgenames)[1]
            inc_matrix[j, idx] = 1
        end
        for i in 1:length(unit.outlets)
            streamname = unit.outlets[i]
            idx = findall(==(streamname), edgenames)[1]
            inc_matrix[j, idx] = -1
        end
    end
    
    return inc_matrix, nodenames, edgenames
end


function signal_matrix(fs::Flowsheet)
    inc_matrix, nodenames, edgenames = incidence_matrix(fs);

    sig_matrix, signal_labels = expand_incidence_matrix(inc_matrix')
    
    return sig_matrix, edgenames, nodenames[signal_labels]
end


function graph_from_incidence_matrix(incmtx, nodenames, edgenames, returnlabels=false)
    numnodes, numedges = size(incmtx)

    graph = DiGraph(numnodes)
    edgetolabel = Dict{Edge, String}()

    for i in 1:numedges
        sources = findall(==(-1), incmtx[:,i])
        destinations = findall(==(1), incmtx[:,i])
        numsources = length(sources)
        numdestinations = length(destinations)

        if numsources == 0 # feed stream
            numnodes += 1
            push!(nodenames, "Source_$(edgenames[i])")
            add_vertex!(graph)
            add_edge!(graph, numnodes, destinations[1])
            edgetolabel[Edge(numnodes, destinations[1])] = edgenames[i]
        elseif numdestinations  == 0 # product stream
            numnodes += 1
            push!(nodenames, "Sink_$(edgenames[i])")
            add_vertex!(graph)
            add_edge!(graph, sources[1], numnodes)
            edgetolabel[Edge(sources[1], numnodes)] = edgenames[i]
        else
            add_edge!(graph, sources[1], destinations[1])
            edgetolabel[Edge(sources[1], destinations[1])] = edgenames[i]
        end
    end

    #TODO Remove the plotting stuff after debugging is over??
    if returnlabels
        return graph, nodenames, edgetolabel
    else
        return graph
    end
end


function expand_incidence_matrix(inc_matrix)

    numrows = size(inc_matrix, 1)
    newmatrix = Int64[]
    unitlabels = []

    firstcol = true
    colnum = 0
    for col in eachcol(inc_matrix)
        colnum += 1

        feeds = findall(==(1), col)
        prods = findall(==(-1), col)
 
        for (outlet, inlet) in Base.Iterators.product(feeds, prods)
            newcol = zeros(Int64, numrows)
            newcol[inlet] = +1
            newcol[outlet] = -1

            if firstcol
                newmatrix = copy(newcol)
                firstcol = false
            else
                newmatrix = hcat(newmatrix, newcol)
            end

            push!(unitlabels, colnum)
        end   
    end     
        
    return newmatrix, unitlabels
end


function graph_reduction(sig_matrix, nodenames, edgenames)
    # Step 1: merge parallel arcs going in the same direction
    # This effectively means removing identical columns from the signal matrix
    # Some extra work to keep track of track the edge names

    unique_cols = unique(eachcol(sig_matrix))
    unique_idxs = [findfirst(x -> x == col, unique_cols) for col in unique_cols]
    sm = reduce(hcat, unique_cols)
    enames = edgenames[unique_idxs]
    nnames = copy(nodenames)

    # Find the nodes with only one preceding node
    numinputs = count(==(1), sm, dims=2)
    todelete = Int64[]
        
    for (node, num_in) in enumerate(numinputs)
        if num_in == 1
            input = findfirst(==(1), sm[node, :])
            source = findfirst(==(-1), sm[:, input])            

            # Merge the nodes
            # Add the downstream node's row in sm to the upstream node's row
            sm[source, :] += sm[node, :]

            # Delete the downstream node's row
            push!(todelete, node)
        end
    end
    sm = sm[setdiff(1:end, todelete), :]
    nnames = nnames[setdiff(1:end, todelete), :]

    # Eliminate nodes with self-loops
    return sm
end







#TODO Remove the plotting stuff after debugging is over??

function generate_graph(fs::Flowsheet, returnlabels=false)
    # Create a directed graph with the unit operations as nodes

    incmtx, nodenames, edgenames = incidence_matrix(fs)
    graph_from_incidence_matrix(incmtx, nodenames, edgenames, returnlabels)
end


function generate_signal_graph(fs::Flowsheet, returnlabels=false)
    # Generate the dual, signal graph where unit ops are edges and streams are nodes

    sigmtx, nodenames, edgenames = signal_matrix(fs)
    graph_from_incidence_matrix(sigmtx, nodenames, edgenames, returnlabels)
end


function plot_graph(fs::Flowsheet)
    graph, nodenames, edgetolabel = generate_graph(fs, true)
    fig, ax, plt = graphplot(graph, nlabels=nodenames, elabels=[edgetolabel[e] for e in edges(graph)])
    hidedecorations!(ax)
    display(fig)
end


function plot_signal_graph(fs::Flowsheet)
    graph, nodenames, edgetolabel = generate_signal_graph(fs, true)
    fig, ax, plt = graphplot(graph, nlabels=nodenames, elabels=[edgetolabel[e] for e in edges(graph)])
    hidedecorations!(ax)
    display(fig)
end









inc_matrix, _, _ = incidence_matrix(fs)
sig_matrix, nodenames, edgenames =  signal_matrix(fs)

inc_matrix
sig_matrix

g = generate_graph(fs)
simplecycles(g)
s = generate_signal_graph(fs)


plot_graph(fs)
plot_signal_graph(fs)