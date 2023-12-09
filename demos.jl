using FlowsheetTools, Statistics


#-Components---------------------------

# We need a ComponentList to hold all the
# components, so we know where to find them later

syscomps = ComponentList()


# You can read them from a folder with saved components
count = readcomponentlist!(syscomps, "components", ["Ethylene", "Ethane", "Hydrogen"])

# Or define them directly

@comp begin
    N --> 2
end "Nitrogen" syscomps

# And then save them to file

writecomponent(joinpath("components/", "Nitrogen.comp"), syscomps["Nitrogen"])


#-Streams------------------------------

#  Like for the components, we need a container for streams
#  so we have something to iterate through later

sysstreams = StreamList()

# You can create the streams directly with instantaneous flows
# This can be in either mass or molar units
# The units are not specified - if you assume the mass flows are
# in kg/h, then the molar equivalent is kmol/hr, but this could
#  as easily be lb/week and lbmole/week.

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

# One stream was specified as mass flows, the other as molar flows,
# but these streams are the same and the missing flows are calculated
# automatically

sysstreams["Test"].moleflows .≈ sysstreams["Product"].moleflows

# As are the atomic flows, in the same units as the molar flows

all(getindex.(values(sysstreams["Test"].atomflows), "C") .== getindex.(values(sysstreams["Product"].atomflows), "C"))
all(getindex.(values(sysstreams["Test"].atomflows), "H") .== getindex.(values(sysstreams["Product"].atomflows), "H"))

#  When we want to deal with streams with multiple historic data points,
#  the best option is to use readstreamhistory to read them from file

sysstreams = StreamList() # Create a new container and dump the previous streams
sysstreams["Feed"] = readstreamhistory(joinpath("streamhistories", "FeedStream.csv"), "Feed", syscomps; ismoleflow=true)
sysstreams["Product"] = readstreamhistory(joinpath("streamhistories", "ProdStream.csv"), "Product", syscomps; ismoleflow=true)

# In the files, we had data for Ethylene, Ethane and Hydrogen, but our list of
# components also includes nitrogen. Zero flows are automatically added for any
# components not in the file, so all streams contain all the components.
# We can also add components after reading the files, but then we need to 
# call refreshcomplist(streamlist) to add the zero flows to all the new components
# to existing streams.

@comp begin
    Ar --> 1
end "Argon" syscomps

refreshcomplist(sysstreams)

sysstreams["Feed"]

# Manipulate some streams

# Multiplication with a scalar
sysstreams["Prod2"] = 2.0*sysstreams["Product"]

# Note the use of .≈ and .* - internally the data are stored in TimeArrays (TimeSeries.jl)
# and only the braodcasted operators are used on TimeArrays. Comparing to TimeArrays also
# returns a TimeArray and we extract the values using values() to get a BitVector
all(values(sysstreams["Prod2"].totalmassflow) .≈ values(2.0 .* sysstreams["Product"].totalmassflow))

# Copy and copy with multiply
copystream!(sysstreams, "Product", "MyStream")
copystream!(sysstreams, "Product", "MyStream2"; factor=2.0) # double the flow!

# Check that they are identical
all(values(sysstreams["MyStream2"].totalmassflow .≈ 2.0 .* sysstreams["MyStream"].totalmassflow))

# Different name, so not identical
sysstreams["Product"] == sysstreams["MyStream"]

# But all the flows are the same
all(getindex.(values(sysstreams["Product"].atomflows), "C") .== getindex.(values(sysstreams["MyStream"].atomflows), "C"))
all(getindex.(values(sysstreams["Product"].atomflows), "H") .== getindex.(values(sysstreams["MyStream"].atomflows), "H"))
all(getindex.(values(sysstreams["Product"].atomflows), "N") .== getindex.(values(sysstreams["MyStream"].atomflows), "N"))

# Do some more things with streams
renamestream!(sysstreams, "MyStream", "Dummy")
deletestream!(sysstreams, "Dummy")


#-UnitOps, Boundaries and KPIs---------

# Test unit ops and mass balance boundaries
# Define some streams
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

# Define some unit ops
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

# Let's split the product a little. We'll need some empty streams
sysstreams["Product1"] = emptystream(sysstreams, "Product1");
sysstreams["Product1a"] = emptystream(sysstreams, "Product1a");
sysstreams["Product1b"] = emptystream(sysstreams, "Product1b");
sysstreams["Product2"] = emptystream(sysstreams, "Product2");
sysstreams["Product3"] = emptystream(sysstreams, "Product3");

# A flow splitter that splits 50% of the product to each of Product1 and Product2
# These streams will have identcal compositions
@unitop begin
    inlets --> ["Product"]
    outlets --> ["Product1", "Product2"]
    calc --> flowsplitter!
    params --> [0.5]
end "ProductSplitter" sysstreams sysunitops
sysunitops["ProductSplitter"]()

# A component splitter that splits Product1 into Product1a and Product1b
# These streams will have different compositions, with 
# the hydrogen split 50:50, 70% of the ethane going to Product1b and the
# remainder of Product1, going to Product1b (the last stream listed)
@unitop begin
    inlets --> ["Product1"]
    outlets --> ["Product1a", "Product1b"]
    calc --> componentplitter!
    params --> Dict([
        "Hydrogen" => Dict(["Product1a" => 0.5]),
        "Ethane" => Dict(["Product1b" => 0.3])
    ])
end "ProductSplitter" sysstreams sysunitops
sysunitops["ProductSplitter"]()

# And then we mix it all again and check that we still have Product
@unitop begin
    inlets --> ["Product1a", "Product1b", "Product2"]
    outlets --> ["Product3"]
    calc --> mixer!
end "Mixer2" sysstreams sysunitops
sysunitops["Mixer2"]()

# Check that the two streams have the same flows
all(values(sysstreams["Product"].massflows .≈ sysstreams["Product3"].massflows))

# Define a mass balance boundary
@boundary begin
    unitops --> ["Mixer", "Reactor"]
end b sysunitops

# Check the closures
b.atomclosures
b.closure
b.total_in.totalmassflow
b.total_out.totalmassflow

# And some defined KPIs
conversion(b, "Ethane")
conversion(b, "Ethylene")
conversion(b, "Hydrogen")
molar_selectivity(b, "Ethylene", "Ethane")



# Let's do the same with some history attached to the streams

sysstreams = StreamList() # Create a new container and dump the previous streams
sysstreams["C2"] = readstreamhistory(joinpath("streamhistories", "C2.csv"), "C2", syscomps; ismoleflow=true)
sysstreams["H2"] = readstreamhistory(joinpath("streamhistories", "Hydrogen.csv"), "H2", syscomps; ismoleflow=true)
sysstreams["Product"] = readstreamhistory(joinpath("streamhistories", "Product.csv"), "Product", syscomps; ismoleflow=true)
sysstreams["Mixed"] = emptystream(sysstreams, "Mixed")

# Define some unit ops
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
end "ProductSplitter" sysstreams sysunitops
sysunitops["ProductSplitter"]()

@unitop begin
    inlets --> ["Product1a", "Product1b", "Product2"]
    outlets --> ["Product3"]
    calc --> mixer!
end "Mixer2" sysstreams sysunitops
sysunitops["Mixer2"]()

# Check that the two streams have the same flows
all(values(sysstreams["Product"].massflows .≈ sysstreams["Product3"].massflows))


# Define a mass balance boundary
@boundary begin
    unitops --> ["Mixer", "Reactor"]
end b sysunitops

# Check the closures
b.atomclosures
b.closure
b.total_in.totalmassflow
b.total_out.totalmassflow

# And some defined KPIs
c1 = conversion(b, "Ethane")
c2 = conversion(b, "Ethylene")
sc2 = molar_selectivity(b, "Ethylene", "Ethane")

mean(values(c1))
mean(values(c2))
mean(values(sc2))


# Lets introduce some errors and check our closure corrections
# copystream!() can take a factor which it multiplies the flows in the source stream with
copystream!(sysstreams, "C2", "eC2", factor = 1.05)
copystream!(sysstreams, "H2", "eH2", factor = 0.95)
copystream!(sysstreams, "Product", "eProduct")
sysstreams["eMixed"] = emptystream(sysstreams, "eMixed") # We'll calculate this stream with the mixer model


# Define the unit ops and boundary
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

# Define a mass balance boundary
@boundary begin
    unitops --> ["eMixer", "eReactor"]
end b sysunitops


# Get the correction factors on the inlets and outlets
corrections = calccorrections(b, "eProduct")
b2 = closemb(b, anchor = "eProduct")
sysunitops["eMixer"]()

mean(values(b.closure))
mean(values(b2.closure))

print(showdata(b2))


# Test Flowsheet

fs = Flowsheet(sysunitops, ["Reactor"], [1])
addunitop!(fs, ["Mixer"])
fs()

generateBFD(fs, "./myflowsheet.svg")
