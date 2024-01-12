using FlowsheetTools, Statistics

syscomps = ComponentList()

count = readcomponentlist!(syscomps, "components", ["Ethylene", "Ethane", "Hydrogen"])

@comp begin
    N --> 2
end "Nitrogen" syscomps

writecomponent(joinpath("components/", "Nitrogen.comp"), syscomps["Nitrogen"])

sysstreams = StreamList()

@stream mass begin
    "Ethylene" --> 2.8053
    "Ethane" --> 27.06192
    "Hydrogen" --> 2.21738
end "Test" syscomps sysstreams

@stream mole begin
    "Ethane" --> 0.9
    "Hydrogen" --> 1.1
    "Ethylene" --> 0.1
end "Product" syscomps sysstreams

sysstreams["Test"].moleflows .≈ sysstreams["Product"].moleflows

sysstreams["Test"] ≈ sysstreams["Product"]

sysstreams["Test"] == sysstreams["Product"]

all(getindex.(values(sysstreams["Test"].atomflows), "C") .== getindex.(values(sysstreams["Product"].atomflows), "C"))

all(getindex.(values(sysstreams["Test"].atomflows), "H") .== getindex.(values(sysstreams["Product"].atomflows), "H"))

sysstreams = StreamList()

sysstreams["Feed"] = readstreamhistory(joinpath("streamhistories", "FeedStream.csv"), "Feed", syscomps; ismoleflow=true)
sysstreams["Product"] = readstreamhistory(joinpath("streamhistories", "ProdStream.csv"), "Product", syscomps; ismoleflow=true)

@comp begin
    Ar --> 1
end "Argon" syscomps

refreshcomplist(sysstreams)

sysstreams["Feed"]

sysstreams["Prod2"] = 2.0*sysstreams["Product"]

sysstreams["Prod2"] == 2.0*sysstreams["Product"]

all(values(sysstreams["Prod2"].totalmassflow) .== values(2.0 .* sysstreams["Product"].totalmassflow))

copystream!(sysstreams, "Product", "MyStream")
copystream!(sysstreams, "Product", "MyStream2"; factor=2.0)
sysstreams["MyStream2"] ≈ 2.0 * sysstreams["MyStream"]

copystream!(sysstreams, "Product", "MyStream")
sysstreams["Product"] == sysstreams["MyStream"]

copystream!(sysstreams, "Product", "MyStream2"; factor=2.0) # double the flow!
2.0 * sysstreams["Product"] == sysstreams["MyStream2"]

all(getindex.(values(sysstreams["Product"].atomflows), "C") .== getindex.(values(sysstreams["MyStream"].atomflows), "C"))

all(getindex.(values(sysstreams["Product"].atomflows), "H") .== getindex.(values(sysstreams["MyStream"].atomflows), "H"))

all(getindex.(values(sysstreams["Product"].atomflows), "N") .== getindex.(values(sysstreams["MyStream"].atomflows), "N"))

renamestream!(sysstreams, "MyStream", "Dummy")
deletestream!(sysstreams, "Dummy")
sysstreams

sysstreams = StreamList()

@stream mole begin
    "Hydrogen" --> 1.1
end "H2" syscomps sysstreams

@stream mole begin
    "Ethylene" --> 0.1
    "Ethane" --> 0.9
end "C2" syscomps sysstreams

sysstreams["Mixed"] = emptystream(sysstreams, "Mixed")

@stream mole begin
    "Ethylene" --> 0.0
    "Ethane" --> 1.0
    "Hydrogen" --> 1.0
end "Product" syscomps sysstreams

sysunitops = UnitOpList()

@unitop begin
    inlets --> ["H2", "C2"]
    outlets --> ["Mixed"]
    calc --> mixer!
end "Mixer" sysstreams sysunitops

sysunitops["Mixer"]()

@unitop begin
    inlets --> ["Mixed"]
    outlets --> ["Product"]
end "Reactor" sysstreams sysunitops

sysstreams["Product1"] = emptystream(sysstreams, "Product1");
sysstreams["Product1a"] = emptystream(sysstreams, "Product1a");
sysstreams["Product1b"] = emptystream(sysstreams, "Product1b");
sysstreams["Product2"] = emptystream(sysstreams, "Product2");
sysstreams["Product3"] = emptystream(sysstreams, "Product3");

@unitop begin
    inlets --> ["Product"]
    outlets --> ["Product1", "Product2"]
    calc --> flowsplitter!
    params --> [0.5]
end "ProductSplitter" sysstreams sysunitops

sysunitops["ProductSplitter"]()

@unitop begin
    inlets --> ["Product1"]
    outlets --> ["Product1a", "Product1b"]
    calc --> componentplitter!
    params --> Dict([
        "Hydrogen" => Dict(["Product1a" => 0.5]),
        "Ethane" => Dict(["Product1b" => 0.3])
    ])
end "ComponentSplitter" sysstreams sysunitops

sysunitops["ComponentSplitter"]()

@unitop begin
    inlets --> ["Product1a", "Product1b", "Product2"]
    outlets --> ["Product3"]
    calc --> mixer!
end "Mixer2" sysstreams sysunitops

sysunitops["Mixer2"]()

sysstreams["Product"] ≈ sysstreams["Product3"]

@boundary begin
    unitops --> ["Mixer", "Reactor"]
end b sysunitops

b.atomclosures

b.closure

b.total_in.totalmassflow

b.total_out.totalmassflow

conversion(b, "Ethane")

(conversion(b, "Ethylene"), conversion(b, "Hydrogen"))

molar_selectivity(b, "Ethylene", "Ethane")

sysstreams = StreamList() # Create a new container and dump the previous streams
sysstreams["C2"] = readstreamhistory(joinpath("streamhistories", "C2.csv"), "C2", syscomps; ismoleflow=true)
sysstreams["H2"] = readstreamhistory(joinpath("streamhistories", "Hydrogen.csv"), "H2", syscomps; ismoleflow=true)
sysstreams["Product"] = readstreamhistory(joinpath("streamhistories", "Product.csv"), "Product", syscomps; ismoleflow=true)
sysstreams["Mixed"] = emptystream(sysstreams, "Mixed");
sysstreams["Product1"] = emptystream(sysstreams, "Product1");
sysstreams["Product1a"] = emptystream(sysstreams, "Product1a");
sysstreams["Product1b"] = emptystream(sysstreams, "Product1b");
sysstreams["Product2"] = emptystream(sysstreams, "Product2");
sysstreams["Product3"] = emptystream(sysstreams, "Product3");

sysunitops = UnitOpList()

@unitop begin
    inlets --> ["H2", "C2"]
    outlets --> ["Mixed"]
    calc --> mixer!
end "Mixer" sysstreams sysunitops
sysunitops["Mixer"]()

@unitop begin
    inlets --> ["Mixed"]
    outlets --> ["Product"]
end "Reactor" sysstreams sysunitops

@unitop begin
    inlets --> ["Product"]
    outlets --> ["Product1", "Product2"]
    calc --> flowsplitter!
    params --> [0.5]
end "ProductSplitter" sysstreams sysunitops
sysunitops["ProductSplitter"]()

@unitop begin
    inlets --> ["Product1"]
    outlets --> ["Product1a", "Product1b"]
    calc --> componentplitter!
    params --> Dict([
        "Hydrogen" => Dict(["Product1a" => 0.5]),
        "Ethane" => Dict(["Product1b" => 0.3])
    ])
end "ComponentSplitter" sysstreams sysunitops
sysunitops["ComponentSplitter"]()

@unitop begin
    inlets --> ["Product1a", "Product1b", "Product2"]
    outlets --> ["Product3"]
    calc --> mixer!
end "Mixer2" sysstreams sysunitops
sysunitops["Mixer2"]()

sysstreams["Product"] ≈ sysstreams["Product3"]

@boundary begin
    unitops --> ["Mixer", "Reactor"]
end b sysunitops

b.atomclosures

b.closure

b.total_in.totalmassflow

b.total_out.totalmassflow

c1 = conversion(b, "Ethane")
c2 = conversion(b, "Ethylene")

sc2 = molar_selectivity(b, "Ethylene", "Ethane")

(mean(values(c1)), mean(values(c2)), mean(values(sc2)))

copystream!(sysstreams, "C2", "eC2", factor = 1.05)
copystream!(sysstreams, "H2", "eH2", factor = 0.95)
copystream!(sysstreams, "Product", "eProduct")
sysstreams["eMixed"] = emptystream(sysstreams, "eMixed") # We'll calculate this stream with the mixer model

@unitop begin
    inlets --> ["eH2", "eC2"]
    outlets --> ["eMixed"]
    calc --> mixer!
end "eMixer" sysstreams sysunitops
sysunitops["eMixer"]()

@unitop begin
    inlets --> ["eMixed"]
    outlets --> ["eProduct"]
end "eReactor" sysstreams sysunitops

@boundary begin
    unitops --> ["eMixer", "eReactor"]
end b sysunitops

corrections = calccorrections(b, "eProduct")

#= `calccorrections` takes a boundary for which to calculate the correction factors, an anchor stream, for which the correction is always 1.0 - no change, and then optional weights for the total mass balance error and the elemental errors.
These latter values default to 1.0 each.

    `function calccorrections(boundary::BalanceBoundary, anchor::String; totalweight=1.0, elementweight=1.0)`
=#

b2 = closemb_simple(b, anchor = "eProduct")

(mean(values(b.closure)), mean(values(b2.closure)))

print(showdata(b2))

fs = Flowsheet(sysunitops, ["Reactor"], [1])
addunitop!(fs, ["Mixer", "ProductSplitter", "ComponentSplitter", "Mixer2"])

fs()

generateBFD(fs, "./myflowsheet.svg", displaybfd=false)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
