using FlowsheetTools

syscomps = ComponentList()
count = readcomponentlist!(syscomps, "components", ["Ethylene", "Ethane", "Hydrogen"])




# The next one should error
# sysstreams["Wrong"] = readstreamhistory(joinpath("streamhistories", "WrongStream.csv"), "Wrong", syscomps; ismoleflow=true)

dummy = emptystream(sysstreams, "Dummy")
dummy = sysstreams["Feed"]
dummies = sysstreams[["Feed", "Product"]]

count = length(sysstreams)

dummy = dummies[1] + dummies[2]
dummy = 2 * dummy
dummy = dummy * 2

shortstreams = StreamList()
shortstreams["Feed"] = readstreamhistory(joinpath("streamhistories", "FeedStreamShort.csv"), "Feed", syscomps; ismoleflow=true)

# Iterate through the CompList - returns pairs of names and Component objects
for (compname, comp) in syscomps
    println(compname)
end

# Iterate through the StreamList - returns pairs of names and Stream objects
for (streamname, stream) in sysstreams
    println(streamname)
end

sysunits = UnitOpList()

@unitop begin
    inlets --> ["Feed"]
    outlets --> ["Product"]
    calc --> mixer!
end "Mixer" sysstreams sysunits

# Iterate through the UnitOpList - returns pairs of names and UnitOp objects
for (unitopname, unitop) in sysunits
    println(unitopname)
end
