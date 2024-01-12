# # FlowsheetTools.jl Demonstration
# FlowsheetTools.jl is a simply library for dealing with flowsheets (components, streams, unitops, boundaries and flowsheets).
# It can be used as a platform for running custom models, for example when fitting kinetic parameters to pilot plant data, where the operating unit is more complicated than a single reactor.
# The primary intended purpose however, was for process analytics - generating KPIs on a flowsheet and reconciling mass balances for generic flowsheets.

using FlowsheetTools, Statistics

# ## Components

# We need a ComponentList to hold all the components, so we know where to find them later

syscomps = ComponentList()

# You can read them from a folder with saved components

count = readcomponentlist!(syscomps, "components", ["Ethylene", "Ethane", "Hydrogen"])

# Or define them directly

@comp begin
    N --> 2
end "Nitrogen" syscomps

# And then save them to file to re-use later.

writecomponent(joinpath("components/", "Nitrogen.comp"), syscomps["Nitrogen"])


# ## Streams

# As for components, we create a container stream list to hold the streams so we have something to iterate through later.

sysstreams = StreamList()

# You can create the streams directly with instantaneous flows.
# This can be in either mass or molar flows.
# The units are not specified - if you assume the mass flows are in kg/h, then the molar equivalent is kmol/hr, but this could as easily be lb/week and lbmole/week.

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

# One stream here was specified as mass flows, the other as molar flows, but there streams are the same and the missing flows (mass/mole) are calculated automatically in the constructor.

# We can quickly check if the molar flows are identical:

sysstreams["Test"].moleflows .≈ sysstreams["Product"].moleflows

# Or, more conveniently, directly with the `≈` or `==` operators.
# Keep in mind that using `==` for floating point values is likely to give `false` when you would expect `true`, so it is recommende to rather use `≈` (`\approx<tab>`)

sysstreams["Test"] ≈ sysstreams["Product"]
#-
sysstreams["Test"] == sysstreams["Product"]

# And, for the skeptical members of the audience, we can also check the atomic flows:

all(getindex.(values(sysstreams["Test"].atomflows), "C") .== getindex.(values(sysstreams["Product"].atomflows), "C"))
#-
all(getindex.(values(sysstreams["Test"].atomflows), "H") .== getindex.(values(sysstreams["Product"].atomflows), "H"))

# When we want to deal with streams with multiple historic data points, we read them from a file.

# First, start with a new, empty stream list:

sysstreams = StreamList() 

# Then we read the streams from file.

sysstreams["Feed"] = readstreamhistory(joinpath("streamhistories", "FeedStream.csv"), "Feed", syscomps; ismoleflow=true)
sysstreams["Product"] = readstreamhistory(joinpath("streamhistories", "ProdStream.csv"), "Product", syscomps; ismoleflow=true)

# In the data files (*.csv), we had columns of data for ethylene, ethane and hydrogen, but or list of components also include nitrogen.
# We automatically set zero flows for amy components not in the file, so all the streams contain all of the components (for our sanity).

# We can still add components to the component list after the streams were created.
# If we do, then we should also call `refreshcomplist(streamlist)` to add zero flows for all of these new components to the existing streams in the stream list.
    
@comp begin
    Ar --> 1
end "Argon" syscomps

refreshcomplist(sysstreams)

sysstreams["Feed"];


# ## What can we do with streams?

# Operations defined on streams include addition and multiplication with a scalar. Addition of streams is effectively a mixer unit.
# Multiplication is used to allow correction factors for mass balance reconciliation.

sysstreams["Prod2"] = 2.0*sysstreams["Product"]

# Check the answer

sysstreams["Prod2"] == 2.0*sysstreams["Product"]

# Alternatively

all(values(sysstreams["Prod2"].totalmassflow) .== values(2.0 .* sysstreams["Product"].totalmassflow))

# Note the use of `.==` and `.*` above. Internally the data are stored in `TimeArrays` from `TimeSeries.jl` and only the broadcasted operators are used on `TimeArray`s.

# Comparison between `TimeArrays` returns a `TimeArray` and we extract the results as an aray using the `values()` function to get a `BitVector`.

# We can also copy streams and copy with a multiplication factor:

copystream!(sysstreams, "Product", "MyStream")
copystream!(sysstreams, "Product", "MyStream2"; factor=2.0)
sysstreams["MyStream2"] ≈ 2.0 * sysstreams["MyStream"]

# Copy and copy with multiply

copystream!(sysstreams, "Product", "MyStream")
sysstreams["Product"] == sysstreams["MyStream"]

copystream!(sysstreams, "Product", "MyStream2"; factor=2.0) # double the flow!
2.0 * sysstreams["Product"] == sysstreams["MyStream2"]

# The streams have different names, but we overload `==` to only check the molar flows of each component, so we get the expected answer.

# Since the atomic flows are automatically calculated, they will also match

all(getindex.(values(sysstreams["Product"].atomflows), "C") .== getindex.(values(sysstreams["MyStream"].atomflows), "C"))
#-
all(getindex.(values(sysstreams["Product"].atomflows), "H") .== getindex.(values(sysstreams["MyStream"].atomflows), "H"))
#-
all(getindex.(values(sysstreams["Product"].atomflows), "N") .== getindex.(values(sysstreams["MyStream"].atomflows), "N"))

# We can also rename or delete streams from the stream list:

renamestream!(sysstreams, "MyStream", "Dummy")
deletestream!(sysstreams, "Dummy")
sysstreams

# # UnitOps, Boundaries and KPIs

# Let's start with an empty stream list

sysstreams = StreamList()

@stream mole begin
    "Hydrogen" --> 1.1
end "H2" syscomps sysstreams

@stream mole begin
    "Ethylene" --> 0.1
    "Ethane" --> 0.9
end "C2" syscomps sysstreams

# We can also add an empty stream, since we don't measure the mixed stream. We'll calculate it with a mixer model later.

sysstreams["Mixed"] = emptystream(sysstreams, "Mixed")

@stream mole begin
    "Ethylene" --> 0.0
    "Ethane" --> 1.0
    "Hydrogen" --> 1.0
end "Product" syscomps sysstreams

# Now we define some unit operations. As with components and streams we need a container to be able to access the streams again later.

sysunitops = UnitOpList()

@unitop begin
    inlets --> ["H2", "C2"]
    outlets --> ["Mixed"]
    calc --> mixer!
end "Mixer" sysstreams sysunitops

# And to execute the unit operation, we simply call it.

sysunitops["Mixer"]()

# This `UnitOp` takes the required inlet and outlet streams, but is also assigned a calculation.
# In this case, it is the predefined `mixer!` function, which is a simple stream mixer.
# This can however be any user-defined function, with the correct form.
# These calculations will supply the contents of the outlet streams based on the inlets streams and supplied model parameters.
# They are only needed if there is no information on the outlet streams.

@unitop begin
    inlets --> ["Mixed"]
    outlets --> ["Product"]
end "Reactor" sysstreams sysunitops

# Our `Reactor` does not have an associated calculation. It is just a node in the flowsheet graph, so we shall need information for all of the inlets and outlets.
#-
# Let's split the product a little. We'll need some empty streams.

sysstreams["Product1"] = emptystream(sysstreams, "Product1");
sysstreams["Product1a"] = emptystream(sysstreams, "Product1a");
sysstreams["Product1b"] = emptystream(sysstreams, "Product1b");
sysstreams["Product2"] = emptystream(sysstreams, "Product2");
sysstreams["Product3"] = emptystream(sysstreams, "Product3");

# A flow splitter that splits 50% of the product to each of Product1 and Product2. These streams will have identcal compositions

@unitop begin
    inlets --> ["Product"]
    outlets --> ["Product1", "Product2"]
    calc --> flowsplitter!
    params --> [0.5]
end "ProductSplitter" sysstreams sysunitops

sysunitops["ProductSplitter"]()

# A component splitter that splits Product1 into Product1a and Product1b.
# These streams will have different compositions, with the hydrogen split 50:50, 70% of the ethane going to Product1b and the remainder of Product1, going to Product1b (the last stream listed).

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

# And then we mix it all again and check that we still have the original Product stream

@unitop begin
    inlets --> ["Product1a", "Product1b", "Product2"]
    outlets --> ["Product3"]
    calc --> mixer!
end "Mixer2" sysstreams sysunitops

sysunitops["Mixer2"]()

# Check that the two streams have the same flows
sysstreams["Product"] ≈ sysstreams["Product3"]

# Mass balances and KPIs are defined on a boundary around a number of unit operations. We therefore define a `Boundary` and list the contained `UnitOp`s

@boundary begin
    unitops --> ["Mixer", "Reactor"]
end b sysunitops

# We can look at total mass and elemental closures, as well as the combined in- and outflows.

b.atomclosures
#-
b.closure
#-
b.total_in.totalmassflow
#-
b.total_out.totalmassflow

# We can also define KPIs on the boundary. Here we use the pre-defined KPIs of `conversion(boundary, component)` and `selectivity(boundary, reactant, product)`

conversion(b, "Ethane")

# Ethane was produced, not consumed, so has a negative value for conversion.

(conversion(b, "Ethylene"), conversion(b, "Hydrogen"))

# We had complete conversion of ethylene and only ~9% of hydrogen, due to the large excess fed.

molar_selectivity(b, "Ethylene", "Ethane")

# All of the reacted ethylene was converted to ethane.

# ### For streams with time series data

# Now we can repeat this for streams with multiple historic data points attached:

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

# Empty the unit operation list as well, so we start fresh.

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

# Check that the two streams have the same flows.

sysstreams["Product"] ≈ sysstreams["Product3"]

# Define the mass balance boundary for closures and KPIs

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

# So, we have average conversions of ethane (-11%, meaning it was produced, not consumed), ethylene (99.9%) and selectivity of ethylene conversion to ethane (~100%) similar to the single data point above.

# ## Mass balance reconciliation

# The mass balance reconciliation algorithm is currently *VERY BASIC*! This will be updated at the first opportunity, but will be invisible to the end-user and will not have major impacts on the user interface unless additional user input is required.

# To demomstrate the use of the reconciliation tool, we repeat the flowsheet above, but introduce some (artificial) flow measurement errors.

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


# We can request the correction factors, without applying them.

corrections = calccorrections(b, "eProduct")

#= `calccorrections` takes a boundary for which to calculate the correction factors, an anchor stream, for which the correction is always 1.0 - no change, and then optional weights for the total mass balance error and the elemental errors.
These latter values default to 1.0 each.

    `function calccorrections(boundary::BalanceBoundary, anchor::String; totalweight=1.0, elementweight=1.0)`
=#

# We can apply the corrections, with `closemb()`, which will either take a `Dict` of correction factors, or calculate them automatically, if not specified.

b2 = closemb_simple(b, anchor = "eProduct")

# Let's compare the raw and reconciled closures:

(mean(values(b.closure)), mean(values(b2.closure)))

# We can also request some information from a bounary. This is given in table form, packed into a string.

#jl showdata(b2)
#nb print(showdata(b2))
print(showdata(b2)) #src


# ## Flowsheets

# Lastly, for convenience, we can creat a `Flowsheet` object, which holds a number of unit operations and an execution order.
# If the flowsheet is then executed, each unit operation is execute in order, as specified.
# Unit operations can be added or deleted with utility functions and the execution order can be modified.

fs = Flowsheet(sysunitops, ["Reactor"], [1])
addunitop!(fs, ["Mixer", "ProductSplitter", "ComponentSplitter", "Mixer2"])

fs()

# Lastly, once a `Flowsheet` object is created, a block flow diagram can also be generated.

#nb generateBFD(fs, "./myflowsheet.svg")
#jl generateBFD(fs, "./myflowsheet.svg", displaybfd=false);
generateBFD(fs, "./myflowsheet.svg") #src