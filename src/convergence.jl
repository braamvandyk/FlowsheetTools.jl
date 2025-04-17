using FlowsheetTools
using Graphs
using GraphMakie, CairoMakie



fs = Flowsheet()
_ = readcomponentlist!(fs, "components", ["Ethylene", "Ethane", "Hydrogen"])

@stream mole begin
    "Hydrogen" --> 1.1
end "H2" fs

@stream mole begin
    "Ethylene" --> 0.1
    "Ethane" --> 0.9
end "C2" fs

@stream mole begin
    "Ethylene" --> 0.0
    "Ethane" --> 1.0
    "Hydrogen" --> 1.0
end "Product" fs

addemptystream!(fs, "Mixed")
addemptystream!(fs, "RXOutlet")
addemptystream!(fs, "Product")
addemptystream!(fs, "Recycle")

@unitop begin
    inlets --> ["H2", "C2", "Recycle"]
    outlets --> ["Mixed"]
    calc --> mixer!
end "Mixer" fs

@unitop begin
    inlets --> ["Mixed"]
    outlets --> ["RXOutlet"]
end "Reactor" fs

@unitop begin
    inlets --> ["RXOutlet"]
    outlets --> ["Product", "Recycle"]
    calc --> flowsplitter!
    params --> [0.8]
end "RecycleSplitter" fs

fs()

# generateBFD(fs, "./teartest.svg")


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

    incidence_matrix = zeros(Int, numnodes, numedges)

    for (j, unitpair) in enumerate(fs.unitops )
        unit = unitpair.second #the UnitOp object
        for i in 1:length(unit.inlets)
            streamname = unit.inlets[i]
            idx = findall(==(streamname), edgenames)[1]
            incidence_matrix[j, idx] = 1
        end
        for i in 1:length(unit.outlets)
            streamname = unit.outlets[i]
            idx = findall(==(streamname), edgenames)[1]
            incidence_matrix[j, idx] = -1
        end
    end
    
    return incidence_matrix, nodenames, edgenames
end


function generate_graph(fs::Flowsheet, returnlabels=false)
    # Create a directed graph with the unit operations as nodes

    incmtx, nodenames, edgenames = incidence_matrix(fs)
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

#TODO Remove the plotting stuff after debugging is over??
function plot_graph(fs::Flowsheet)
    graph, nodenames, edgetolabel = generate_graph(fs, true)
    gplt = graphplot(graph, nlabels=nodenames, elabels=[edgetolabel[e] for e in edges(graph)])
    display(gplt)
end

inc_matrix, _, _ = incidence_matrix(fs);
inc_matrix

plot_graph(fs)
g = generate_graph(fs)
simplecycles(g)