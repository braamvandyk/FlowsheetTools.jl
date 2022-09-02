using FlowsheetTools


#-Components---------------------------


syscomps = ComponentList()

@comp begin
    H --> 2
end "Hydrogen" syscomps

@comp begin
    C --> 2
    H --> 6
end "Ethane" syscomps

@comp begin
    C --> 2
    H --> 4
end "Ethylene" syscomps


#-Streams------------------------------


sysstreams = StreamList()

@stream mass begin
    "Ethylene" --> 2.8053
    "Ethane" --> 27.06192
    "Hydrogen" --> 2.21738
end "Test" syscomps sysstreams

@stream mole begin 
    "Ethylene" --> 0.1
    "Ethane" --> 0.9
    "Hydrogen" --> 1.1
end "Product" syscomps sysstreams

sysstreams["Test"].moleflows ≈ sysstreams["Product"].moleflows

# Manipulate some streams

# Copy and copy with multiply
copystream!(sysstreams, "Product", "mystream")
copystream!(sysstreams, "Product", "mystream2"; factor=2.0) # double the flow!

sysstreams["mystream2"].totalmassflow ≈ 2.0 * sysstreams["mystream"].totalmassflow

# Different name, so not identical
sysstreams["Product"] == sysstreams["mystream"]

# But all the flows are the same
sysstreams["Product"].atomflows == sysstreams["mystream"].atomflows

# Do some more things with streams
renamestream!(sysstreams, "mystream", "dummy")
deletestream!(sysstreams, "dummy")

# Multiplication with a scalar
sysstreams["Prod2"] = 2.0*sysstreams["Product"]
sysstreams["Prod2"].totalmassflow ≈ 2.0 * sysstreams["Product"].totalmassflow


#-UnitOps, Boundaries and KPIs---------


# Test unit ops and mass balance boundaries
# Define some streams
sysstreams = StreamList()

@stream mole begin
    "Ethylene" --> 1.0
    "Hydrogen" --> 2.0
end "Feed" syscomps sysstreams

@stream mole begin
    "Ethylene" --> 0.1
    "Ethane" --> 0.9
    "Hydrogen" --> 1.1
end "Product" syscomps sysstreams


@stream mole begin
    "Hydrogen" --> 1.1
end "H2" syscomps sysstreams

@stream mole begin
    "Ethylene" --> 0.1
    "Ethane" --> 0.9
end "C2" syscomps sysstreams

@stream mole begin
    "Hydrogen" --> 0.0
end "Mixed" syscomps sysstreams

# Define some unit ops
sysunitops = UnitOpList()

@unitop begin
    inlets --> ["Feed"]
    outlets --> ["Product"]
end "Reactor" sysstreams sysunitops

@unitop begin
    inlets --> ["Product"]
    outlets --> ["C2", "H2"]
end "Membrane" sysstreams sysunitops

@unitop begin
    inlets --> ["H2", "C2"]
    outlets --> ["Mixed"]
    calc --> mixer!
end "Mixer" sysstreams sysunitops
sysunitops["Mixer"]()

# Define a mass balance boundary
@boundary begin
    unitops --> ["Reactor", "Membrane"]
end b sysunitops

# Check the closures
b.atomclosures
b.closure
b.total_in.totalmassflow
b.total_out.totalmassflow

# And some defined KPIs
conversion(b, "Ethane")
conversion(b, "Ethylene")
selectivity(b, "Ethylene", "Ethane")

# Now test the reconciliations
# Change some flows to simulate measurement errors
copystream!(sysstreams, "Feed", "Feed2"; factor=0.95)
copystream!(sysstreams, "Product", "Prod2"; factor=1.01)

# Define the unit ops and boundary
@unitop begin
    inlets --> ["Feed2"]
    outlets --> ["Prod2"]
end "Reactor2" sysstreams sysunitops

@unitop begin
    inlets --> ["Prod2"]
    outlets --> ["C2", "H2"]
end "Membrane2" sysstreams sysunitops

@boundary begin
    unitops --> ["Reactor2", "Membrane2"]
end b2 sysunitops

# Get the correction factors on the inlets and outlets
corrections = calccorrections(b2)
b2 = closemb(b2)
b2

# Test the pretty printing
print(syscomps["Ethane"])
print(sysstreams["Feed"])
print(sysunitops["Membrane"])
print(b)

# Read in a list of components
syscomps2 = ComponentList()
count = readcomponentlist!(syscomps2, "components", ["Ethylene", "Ethane", "Hydrogen"])
count == 3


#-Streams with historical data---------


histstreams = StreamHistoryList()
# Read in some stream histories
histstreams["Feed"] = readstreamhistory(joinpath("streamhistories", "FeedStream.csv"), "Feed", syscomps; ismoleflow=true)
histstreams["Product"] = readstreamhistory(joinpath("streamhistories", "ProdStream.csv"), "Product", syscomps; ismoleflow=true)
histstreams
# And mix them
histstreams["Comb"] = histstreams["Feed"] + histstreams["Product"]

feeddata = showdata(histstreams["Feed"]);
println(feeddata)
# And scale the result
histstreams["Comb2"] = 2.0 * histstreams["Comb"]


#-UnitOps etc with historical data-----

# Test UnitOpHistory
histops = UnitOpHistoryList()

@unitophist begin
    inlets --> ["Feed", "Product"]
    outlets --> ["Comb2"]
    calc --> mixer!
end "Mixer" histstreams histops
histops["Mixer"]()

histstreams["Comb"].massflows == histstreams["Comb2"].massflows
histstreams["Comb"].comps == histstreams["Comb2"].comps
histstreams["Comb"].timestamps == histstreams["Comb2"].timestamps


# Test BalanceBoundaryHistory
@unitophist begin
    inlets --> ["Feed"]
    outlets --> ["Product"]
end "RX101" histstreams histops

@boundaryhist begin
    unitops --> ["RX101"]
end bh histops

corrections = calccorrections(bh)
bh = closemb(bh)
conversion(bh, "Ethylene")
selectivity(bh, "Ethylene", "Ethane")

print(showdata(bh))


# Test Flowsheet

sysunitops = UnitOpList()

@unitop begin
    inlets --> ["Feed"]
    outlets --> ["Product"]
end "Reactor" sysstreams sysunitops

@unitop begin
    inlets --> ["Product"]
    outlets --> ["C2", "H2"]
end "Membrane" sysstreams sysunitops

@unitop begin
    inlets --> ["H2", "C2"]
    outlets --> ["Mixed"]
    calc --> mixer!
end "Mixer" sysstreams sysunitops


fs = Flowsheet([sysunitops["Reactor"]], [1])
addunitop!(fs, [sysunitops["Membrane"], sysunitops["Mixer"]])
fs()

sysstreams["Dummy"] = sysstreams["H2"] + sysstreams["C2"]
sysstreams["Dummy"].massflows == sysstreams["Mixed"].massflows
sysstreams["Dummy"].atomflows == sysstreams["Mixed"].atomflows