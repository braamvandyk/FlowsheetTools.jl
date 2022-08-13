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

@stream "mass" begin
    "Ethylene" --> 2.8053
    "Ethane" --> 27.06192
    "Hydrogen" --> 2.21738
end syscomps "Test" sysstreams

@stream "mole" begin 
    "Ethylene" --> 0.1
    "Ethane" --> 0.9
    "Hydrogen" --> 1.1
end syscomps "Product" sysstreams

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

@stream "mole" begin
    "Ethylene" --> 1.0
    "Hydrogen" --> 2.0
end syscomps "Feed" sysstreams

@stream "mole" begin
    "Ethylene" --> 0.1
    "Ethane" --> 0.9
    "Hydrogen" --> 1.1
end syscomps "Product" sysstreams


@stream "mole" begin
    "Hydrogen" --> 1.1
end syscomps "H2" sysstreams

@stream "mole" begin
    "Ethylene" --> 0.1
    "Ethane" --> 0.9
end syscomps "C2" sysstreams

@stream "mole" begin
    "Hydrogen" --> 0.0
end syscomps "Mixed" sysstreams

# Define some unit ops
sysunitops = UnitOpList()

sysunitops["Reactor"] = UnitOp("RX101", sysstreams, ["Feed"], ["Product"])
sysunitops["Membrane"]  = UnitOp("MX101", sysstreams, ["Product"], ["C2", "H2"])
sysunitops["Mixer"] = UnitOp("MX101", sysstreams, ["C2", "H2"], ["Mixed"], mixer!)
sysunitops["Mixer"]()

# Define a mass balance boundary
b = BalanceBoundary(sysunitops, ["Reactor", "Membrane"])

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
sysunitops["Reactor2"] = UnitOp("RX101", sysstreams, ["Feed2"], ["Prod2"])
sysunitops["Membrane2"] = UnitOp("MX101", sysstreams, ["Prod2"], ["C2", "H2"])
b2 = BalanceBoundary(sysunitops, ["Reactor2", "Membrane2"])

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
histstreams["Feed"] = readstreamhistory(joinpath("streamhistories", "FeedStream.csv"), "Feed", syscomps)
histstreams["Product"] = readstreamhistory(joinpath("streamhistories", "ProdStream.csv"), "Product", syscomps)

# And mix them
histstreams["Comb"] = histstreams["Feed"] + histstreams["Product"]

# And scale the result
histstreams["Comb"] = 2.0 * histstreams["Comb"]

#-UnitOps etc with historical data-----

# Test UnitOpHistory
histops = UnitOpHistoryList()
histops["RX101"] = UnitOpHistory("RX101", histstreams, ["Feed"], ["Product"])

# Test BalanceBoundaryHistory
bh = BalanceBoundaryHistory(histops, ["RX101"])
corrections = calccorrections(bh)
bh = closemb(bh)
bh