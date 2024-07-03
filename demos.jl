# # FlowsheetTools.jl Demonstration
# FlowsheetTools.jl is a library for dealing with flowsheets (components, streams, unit operations, mass balance boundaries, and flowsheets).
# It can be used as a platform for running custom models, for example when fitting kinetic parameters to pilot plant data, where the operating unit is more complicated than a single reactor. The primary intended purpose however, was for process analytics - generating KPIs on a flowsheet and reconciling mass balances for generic flowsheets.

# For more convenient analysis of flowsheets with missing measurements, a few utility unit operations are provided: a mixer, a flow splitter, a component splitter (to emulate a separation process with split factors) and a stoichiometric reactor block with specified conversions. Custom unit operations can also be defined, by simply providing a function that calculates the outlet streams from the inlets and an optional list of parameters.

# The intention was not to build a full-on process simulator, but the custom reactor blocks etc can be easily added, when needed.

# Let's have a look at how to use the library.

using FlowsheetTools, Statistics

# # The Flowsheet

# The `Flowsheet` is the central object in FlowsheetTools.jl. It contains a list of unit operations, a list of streams, a list of mass balance boundaries, and a list of components. You need only create the flowsheet. It will manange the other components automatically.
# Let's create an empty flowsheet and then build it out.

fs = Flowsheet()

# As you can see, the flowsheet is empty. It has no unit operations, streams, mass balance boundaries, or components yet. You will also see an empty execution order. It is possible to ask the flowsheet to execute all of the unitops it comtains. The execution order sets the order in which they execute.
# Remember, that this is not a process simulator. There is no flowsheet convergence algorithm. The unit operations are of the types needed to allow missing information to be calculated to close a mass balance.
# These include mixers, various splitters and a simple stoichiometric reactor.

# ## Components
# The most basic building block we need is a set of components. A component in FlowsheetTools.jl is a fairly simple object. It has a name, so we can refer to it, and contains the list of atoms and the number of each that make up the component. The molar mass of the component is automatically calculated and also stored.

# We store all the components in a `ComponentList`, so we have a container to find them in, when needed.
# You don't need to directly create and manage a `ComponentList`. The `Flowsheet` manages the `ComponentList` for you.
# `ComponentList` is a wrapper around a Dict{String, Component} and can be in the same way as such any `Dict` in Julia, using the component names to index. Now that we have a container to put them into, we can add some components.

# Let's look at the component list we have in our flowsheet:

fs.comps

# Not much to see. Just an empty list. Let's add some components.

# Since it is likely that we shall re-use components, they can be stored in files, so we don't need to define them every time. Let's read in some components created earlier and stored in the sub-folder `components` under the active folder:

count = readcomponentlist!(fs, "components", ["Ethylene", "Ethane", "Hydrogen"])

# The function `readcomponents` returns the number of components read in - 3 in this case. We specified the names of the components to read in from the folder. There can be any number of files stored there. We also supplied the flowsheet (`fs`) as the first argument.
# Since reading in components modifies (mutates) the flowsheet, the function name is ended in an exclamation mark, as is the convention in Julia.

fs.comps

# Now we have three components in our flowsheet.
# To access a component in the component list, we index using the name of the component.

fs.comps["Ethylene"]

# If we need to define new components, ther are two ways to go about this.

# The first is by using the @comp macro. It takes a list of atoms and a number of each, separated by a -->. We also suply a name for the component, so we can find it in the `ComponentList` it will be stored in, and the name of the `Flowsheet` to which this `ComponentList` belongs.

@comp begin
    N --> 2
end "Nitrogen" fs

# And we can check that is was added to our flowsheet:

fs.comps

# The second way is to create the components by calling the constructor directly. This is most useful when creating lists of components, such a homologous series. For example, we could create the n-paraffins from C1 to C10 as follows: 

cl = ComponentList() # Create a new component list that is NOT in our flowsheet. We are throwing this away when we are done.
for n in 1:10
    if n == 1
        name = "CH4"
    else
        name = "C$(n)H$(2n+2)"
    end
    cl[name] = Component(name, ["C", "H"], [n, 2n+2])
end
cl
#-
cl["C10H22"]

# Here we created a "dummy" component list, so as to not add these components to our flowsheet. If we wanted to add them to the flowsheet, the code would look like this:

for n in 1:10
    if n == 1
        name = "CH4"
    else
        name = "C$(n)H$(2n+2)"
    end
    fs.comps[name] = Component(name, ["C", "H"], [n, 2n+2])
end

# As you can see, the internal `ComponentList` in our flowsheet is called `fs.comps`, where `fs` is of course the name of our flowsheet.

# We can save the components to file to re-use later. `writecomponents` takes the path to write to, and the specific component to write. It returns the number of bytes written.

writecomponent(joinpath("components/", "Nitrogen.comp"), fs.comps["Nitrogen"])


# ## Streams

# Now that we have components, we can create streams for our process. Each stream contains a list of components and their flowrates. You can specify either the mass or moalr flows when creating the stream and the other will be automatically calculated. The constructor will also calculate the flowrates for each type of atom in the stream.

# There are of course two ways in which we would use streams. Either with a single, current value for the flowrate and composition, or with a set of historical values of these. The former is useful for simulations, while the latter is useful for analysis. In either case, the flows are stored in `TimeArrays` from `TimeSeries.jl`. In cases where we only have a single flowrate, this is a simply `TimeArray` of length 1, with an arbitrary (zero) timestamp asigned to the value.

# As was the case with components, we need a container (a stream list) to hold the streams so we have something to iterate through later. Just like `ComponentList`, `StreamList` is a wrapper around a `Dict{String, Stream}`.

sysstreams = StreamList()

# We can create the streams directly with instantaneous flows. This can be in either mass or molar flows. While we could call the constructor directly, using the `@stream` macro is more convenient.

# The first parameter indicates whether the flows are mass or molar flows. Similarly to what we did for components with `@comp`, we then provide a list of components and their flowrates, separated with a -->. We also supply a name for the stream, against which it is stored in the `StreamList`, the component list in which the components are defined, and lastly the `StreamList` to which to add the new stream.

# The units are not specified - if you assume the mass flows are in kg/h, then the molar equivalent is kmol/hr, but this could as easily be lb/week and lbmole/week.

# Here we specify two stream, of identical composition and flows, but by specifying mass flows for the first and molar flows for the second.

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

# These comparisons also showed that we can access the flows using the `massflows`, `molarflows` and `atomflows` properties.

sysstreams["Product"].massflows
#-
sysstreams["Product"].moleflows
#-
sysstreams["Product"].atomflows

# When we want to deal with streams with multiple historic data points, to analyse plant data, we can use the `readstreamhistory` function, to read the stream from a file.

# First, let's start with a new, empty stream list, to get rid of the streams we have created earlier:

sysstreams = StreamList() 

# Then we in read the streams. We need to specify the folderpath to the CSV files, the name of the stream, the component list in which the components are defined and the stream list to which to add the new stream. We also specify whether the flows in the CSV files are mass or molar flows. `ismoleflow` has a default value of false, so need not be specified for mass flows.

sysstreams["Feed"] = readstreamhistory(joinpath("streamhistories", "FeedStream.csv"), "Feed", syscomps; ismoleflow=true)
sysstreams["Product"] = readstreamhistory(joinpath("streamhistories", "ProdStream.csv"), "Product", syscomps; ismoleflow=true)

# In the data files (*.csv), we had columns of data for ethylene, ethane and hydrogen, but or list of components also include nitrogen.
# We automatically set zero flows for amy components not in the file, so all the streams contain all of the components (for our sanity).

# In the data files (*.csv), we had columns of data for ethylene, ethane and hydrogen, but our list of components also include nitrogen. We automatically set zero flows for amy components not in the file, so all the streams contain all of the components (for our sanity).

# We can still add components to the component list after the streams were created. If we do, then we should also call `refreshcomplist(streamlist)` to add zero flows for all of these new components to the existing streams in the stream list.
    
@comp begin
    Ar --> 1
end "Argon" syscomps

refreshcomplist(sysstreams)

sysstreams["Feed"]

# ## What can we do with streams?

# Operations defined on streams include adding streams together and multiplying a stream with a scalar value. Addition of streams is effectively a mixer unit.
# Multiplication is used to allow correction factors for mass balance reconciliation.

sysstreams["Prod2"] = 2.0*sysstreams["Product"]

# Let's check the answer:

sysstreams["Prod2"] == 2.0*sysstreams["Product"]

# Alternatively,

all(values(sysstreams["Prod2"].totalmassflow) .== values(2.0 .* sysstreams["Product"].totalmassflow))

# Note the use of `.==` and `.*` above. Internally the data are stored in `TimeArrays` from `TimeSeries.jl` and only the broadcasted operators are used on `TimeArray`s.

# Comparison between `TimeArrays` returns a `TimeArray` with the comparison for each timestamp and we extract the results as an aray using the `values()` function to get a `BitVector`.

# We can also copy streams and combine the two streams by copy with a scalar multiplication in a single call:

copystream!(sysstreams, "Product", "MyStream")
copystream!(sysstreams, "Product", "MyStream2"; factor=2.0) # double the flow!
sysstreams["MyStream2"] ≈ 2.0 * sysstreams["MyStream"]

# The streams have different names, but we overload `==` to only check the molar flows of each component, so we get the expected answer.

# Since the atomic flows are automatically calculated, they will also match

all(getindex.(values(sysstreams["Product"].atomflows), "C") .== getindex.(values(sysstreams["MyStream"].atomflows), "C"))
#-
all(getindex.(values(sysstreams["Product"].atomflows), "H") .== getindex.(values(sysstreams["MyStream"].atomflows), "H"))
#-
all(getindex.(values(sysstreams["Product"].atomflows), "N") .== getindex.(values(sysstreams["MyStream"].atomflows), "N"))

# We can also rename or delete streams from the stream list:

renamestream!(sysstreams, "MyStream", "Dummy")
sysstreams
#- 
deletestream!(sysstreams, "Dummy")
sysstreams

# # UnitOps, Boundaries and KPIs

# Let's start with an empty stream list again

sysstreams = StreamList()

# We'll add some instantaneous flow streams.

@stream mole begin
    "Hydrogen" --> 1.1
end "H2" syscomps sysstreams

@stream mole begin
    "Ethylene" --> 0.1
    "Ethane" --> 0.9
end "C2" syscomps sysstreams

@stream mole begin
    "Ethylene" --> 0.0
    "Ethane" --> 1.0
    "Hydrogen" --> 1.0
end "Product" syscomps sysstreams

# We can also add an empty stream as a placeholder. We'll calculate it with a mixer model later.

sysstreams["Mixed"] = emptystream(sysstreams, "Mixed")

# Now we define some unit operations. As with components and streams we create a container to be able to conveniently access them again later. A `UnitOpList` works the same way as a `ComponentList` or `StreamList` - it is a wrapper around a `Dict{String, UnitOp}`.

sysunitops = UnitOpList()

@unitop begin
    inlets --> ["H2", "C2"]
    outlets --> ["Mixed"]
    calc --> mixer!
end "Mixer" sysstreams sysunitops

# The `@unitop` macro creates a `UnitOp` object and adds it to the `UnitOpList`. We can then refer to it by its name, like we do with `Component` and `Stream` objects.
# The macro takes an array of `Stream` names for inlets, and another for outlets. The `calc` field is optional. If we are only calculating KPIs for our process or reconciling mass balances, unit operations do not need to do calculations. They only serve as nodes where streams are connected.
# If we want to do calculations, like for mixers and splitters, we need to specifiy the name of the function to call in the `calc` field.
# The function takes the form:
#     function functionname!(streamlist::StreamList, outlets::Vector{String}, inlets::Vector{String}, params)
# As per Julia convention, we add a `!` at the end of the function name, which means the function modifies some of the variables passed - the outlet(s).
# There can be multiple inlets and outlets, but even single inlets and outlets must be passed in an array.
# All streams passed to a unit operation must be in the stream list. And all streams in a stream list must use the same component list. This keeps things consistent.
# The `params` field does not have a defined type and can be anything needed for the calculation. We'll look at this in more detail later.


# To execute the unit operation, we simply call it.

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

# Let's also add a dummy stream with a fixed composition. We don't need it here, but it is easy to do
sysstreams["Dummy"] = fixedstream(sysstreams, "Dummy", [10.0, 0.0, 0.0, 0.0, 0.1])

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

sysboundaries = BoundaryList()

@boundary begin
    unitops --> ["Mixer", "Reactor"]
end "B1" sysunitops sysboundaries

# We can look at total mass and elemental closures, as well as the combined in- and outflows.

sysboundaries["B1"].atomclosures
#-
sysboundaries["B1"].closure
#-
sysboundaries["B1"].total_in.totalmassflow
#-
sysboundaries["B1"].total_out.totalmassflow

# We can also define KPIs on the boundary. Here we use the pre-defined KPIs of `conversion(boundary, component)` and `selectivity(boundary, reactant, product)`

conversion(sysboundaries["B1"], "Ethane")

# Ethane was produced, not consumed, so has a negative value for conversion.

(conversion(sysboundaries["B1"], "Ethylene"), conversion(sysboundaries["B1"], "Hydrogen"))

# We had complete conversion of ethylene and only ~9% of hydrogen, due to the large excess fed.

molar_selectivity(sysboundaries["B1"], "Ethylene", "Ethane")

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

# And just to show what happens when we create a fixed composition stream when the streams in the StreamList have time series data
sysstreams["Dummy"] = fixedstream(sysstreams, "Dummy", [10.0, 0.0, 0.0, 0.0, 0.1])

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

sysboundaries = BoundaryList()
@boundary begin
    unitops --> ["Mixer", "Reactor"]
end "B1" sysunitops sysboundaries

sysboundaries["B1"].atomclosures

sysboundaries["B1"].closure

sysboundaries["B1"].total_in.totalmassflow

sysboundaries["B1"].total_out.totalmassflow

c1 = conversion(sysboundaries["B1"], "Ethane")
c2 = conversion(sysboundaries["B1"], "Ethylene")

sc2 = molar_selectivity(sysboundaries["B1"], "Ethylene", "Ethane")

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
end "B2" sysunitops sysboundaries


# We can request the correction factors, without applying them.
corrections = calccorrections(sysboundaries)


# We can also pass through a custom additional error function to minimize.
# We have a custom function that receives a dictionary of the correction factors and can do any arbitrary calculation.
# The retured value is added to the weighted square errors of the total mass and element balances.
# You can, for example, add distance from equilibrium for a reaction to the reconciliation.

myerr(factorDict) = 100.0 * sum(abs2, 1.0 .- values(factorDict))
corrections = calccorrections(sysboundaries, customerror=myerr)

# The previous example uses a constant weight for all elements, but we can specifiy individual weights as well.

weights = Dict(["H" => 1.0, "C" => 1.5, "O" => 1.0, "Ar" => 0.0, "N" => 0.0])
corrections2 = calccorrections(sysboundaries, customerror=myerr, setelements=true, elementweights=weights)

# `calccorrections` takes a boundary for which to calculate the correction factors, and then optional weights for the total mass balance error and the elemental errors.
# These latter values default to 1.0 each. It uses [Ridge Regression](https://www.ibm.com/topics/ridge-regression?_sm_au_=iF5HM658VZnjn6Sr0GLqHKHB1jtC6) with a default λ = 0.1
# 
#   function calccorrections_anchor(boundary::BalanceBoundary, anchor::String; totalweight=1.0, elementweight = 1.0, 
#         setelements = false, elementweights::Dict{String, Float64} = Dict{String, Float64}())

# An alternative option is to use an anchor stream, rather than Ridge Regression. The correction factor for the anchor stream will be 1.0, i.e. there is no adjustment.

# corrections_a = calccorrections_anchor(sysboundaries["B2"], "eProduct")
# corrections_a = calccorrections_anchor(sysboundaries["B2"], "eProduct", myerr)

# Or again, with inidividual element weights

# corrections_a = calccorrections_anchor(sysboundaries["B2"], "eProduct", setelements=true, elementweights=weights)

# `calccorrections_anchor` takes a boundary for which to calculate the correction factors, an anchor stream, for which the correction is always 1.0 - no change, and then optional weights for the total mass balance error and the elemental errors.
# These latter values default to 1.0 each.
# 
#   `function calccorrections(boundary::BalanceBoundary, anchor::String; totalweight=1.0, elementweight=1.0)`
#

# We can apply the corrections, with `closemb()`, which will take a `Dict` of correction factors.

# b2 = closemb(sysboundaries["B2"], corrections)
# b3 = closemb(sysboundaries["B2"], corrections_a)


# Let's compare the raw and reconciled closures:

# (mean(values(b.closure)), mean(values(b2.closure)), mean(values(b3.closure)))

# We can now write corrected streams back to file

# writestreamhistory(sysstreams["C2"], "corrected.csv")

# We can also request some information from a bounary. This is given in table form, packed into a string.

#jl showdata(b2)
# #  nb print(showdata(b2))
# print(showdata(b2)) #src


# ## Flowsheets

# Lastly, for convenience, we can creat a `Flowsheet` object, which holds a number of unit operations and an execution order.
# If the flowsheet is then executed, each unit operation is execute in order, as specified.
# Unit operations can be added or deleted with utility functions and the execution order can be modified.

# fs = Flowsheet(sysunitops, ["Reactor"], [1])
# addunitop!(fs, ["Mixer", "ProductSplitter", "ComponentSplitter", "Mixer2"])

#jl fs();
# #nb fs()
# fs() #src

# Lastly, once a `Flowsheet` object is created, a block flow diagram can also be generated.

# #nb generateBFD(fs, "./myflowsheet.svg")
#jl generateBFD(fs, "./myflowsheet.svg", displaybfd=false);
# generateBFD(fs, "./myflowsheet.svg") #src