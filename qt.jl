using FlowsheetTools

fs = Flowsheet();

count = readcomponentlist!(fs, "components", ["Ethylene", "Ethane", "Hydrogen", "Nitrogen", "Argon"]);


readstreamhistory!(fs, "C2", joinpath("streamhistories", "C2.csv"); ismoleflow=true);
readstreamhistory!(fs, "H2", joinpath("streamhistories", "Hydrogen.csv"); ismoleflow=true);
readstreamhistory!(fs, "Product", joinpath("streamhistories", "Product.csv"); ismoleflow=true);
addemptystream!(fs, "Mixed");
addemptystream!(fs, "Product1");
addemptystream!(fs, "Product1a");
addemptystream!(fs, "Product1b");
addemptystream!(fs, "Product2");
addemptystream!(fs, "Product3");
addemptystream!(fs, "Product4");

addfixedstream!(fs, "Dummy", [10.0, 0.0, 0.0, 0.0, 0.1]);

@unitop begin
    inlets --> ["H2", "C2"]
    outlets --> ["Mixed"]
    calc --> mixer!
end "Mixer" fs;
fs.unitops["Mixer"]();

@unitop begin
    inlets --> ["Mixed"]
    outlets --> ["Product"]
end "Reactor" fs;

@unitop begin
    inlets --> ["Product"]
    outlets --> ["Product1", "Product2"]
    calc --> flowsplitter!
    params --> [0.5]
end "ProductSplitter" fs;
fs.unitops["ProductSplitter"]();

@unitop begin
    inlets --> ["Product1"]
    outlets --> ["Product1a", "Product1b"]
    calc --> componentplitter!
    params --> Dict([
        "Hydrogen" => Dict(["Product1a" => 0.5]),
        "Ethane" => Dict(["Product1b" => 0.3])
    ])
end "ComponentSplitter" fs;
fs.unitops["ComponentSplitter"]();

@unitop begin
    inlets --> ["Product1a", "Product1b", "Product2"]
    outlets --> ["Product3"]
    calc --> mixer!
end "Mixer2" fs;
fs.unitops["Mixer2"]();

fs.streams["Product"] ≈ fs.streams["Product3"];

EthyleneHydrogenation(frac) = Reaction(fs, ["Ethylene", "Hydrogen"], ["Ethane"], [1.0, 1.0], [1.0], "Ethylene", frac);

@unitop begin
    inlets --> ["Product3"]
    outlets --> ["Product4"]
    calc --> stoichiometric_reactor!
    params --> [EthyleneHydrogenation(0.5)]
end "Reactor2" fs;
fs.unitops["Reactor2"]();


# generateBFD(fs, "./myflowsheet.svg")

fs.streams["Mixed"] = 1.1*fs.streams["Mixed"];



@boundary begin
    unitops --> ["Mixer"]
end "B1" fs;

@boundary begin
    unitops --> ["Reactor", "ProductSplitter"]
end "B2" fs;

fs.boundaries["B1"];

corrections = calccorrections(fs, λ = 0.0, anchor = "H2");
closemb!(fs, corrections);

fs.boundaries["B1"];


