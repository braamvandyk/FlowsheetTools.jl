using FlowsheetTools

syscomps = ComponentList()
count = readcomponentlist!(syscomps, "components", ["Ethylene", "Ethane", "Hydrogen", "Nitrogen", "Argon"])


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

sysstreams["Dummy"] = fixedstream(sysstreams, "Dummy", [10.0, 0.0, 0.0, 0.0, 0.1])

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


# fs = Flowsheet(sysunitops, ["Mixer", "Reactor", "ProductSplitter", "ComponentSplitter", "Mixer2"], [1, 2, 3, 4, 5])
# generateBFD(fs, "./myflowsheet.svg")

sysstreams["Mixed"] = 1.1*sysstreams["Mixed"]




sysboundaries = BoundaryList()

@boundary begin
    unitops --> ["Mixer"]
end "B1" sysunitops sysboundaries

@boundary begin
    unitops --> ["Reactor", "ProductSplitter"]
end "B2" sysunitops sysboundaries

sysboundaries["B1"]
sysboundaries["B2"]

corrections = calccorrections(sysboundaries; λ = 0.0, anchor = "H2")
# closemb!(sysboundaries, corrections)

sysboundaries["B1"]
sysboundaries["B2"]