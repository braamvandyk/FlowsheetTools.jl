function readconfig(filename)
    filedata = read(filename, String)
    config = TOML.parse(filedata)
    return config
end

function generate_flowsheet(config)
    fs = Flowsheet()

    # Read components
    compsfolder = config["comps"]["readfromfile"]["folder"]
    compnames = config["comps"]["readfromfile"]["names"]
    if length(compnames) > 0
        count = readcomponentlist!(fs, compsfolder, compnames)
    end
    
    # Define new components
    compsnames = config["comps"]["define"]["names"]
    for name in compsnames
        atoms = config["comps"]["define"][name]["atoms"]
        counts = config["comps"]["define"][name]["counts"]
        fs.comps[name] = Component(name, atoms, counts)
    end

    # Read stream histories
    streamfolder = config["streams"]["readfromfile"]["folder"]
    streamnames = config["streams"]["readfromfile"]["names"]
    streamfiles = config["streams"]["readfromfile"]["filenames"]
    ismoleflows = config["streams"]["readfromfile"]["ismoleflow"]
    if length(streamnames) > 0
        count = readstreamlist!(fs, streamfolder, streamnames, streamfiles, ismoleflows)
    end

    # Create empty streams
    emptystreamnames = config["streams"]["emptystreams"]["names"]
    for name in emptystreamnames
        addemptystream!(fs, name)
    end
    # return fs
    
    # Create fixed composition streams
    fixedstreamnames = config["streams"]["fixedstreams"]["names"]
    for name in fixedstreamnames
        flows = config["streams"]["fixedstreams"][name]["flows"]
        ismoleflow = config["streams"]["fixedstreams"][name]["ismoleflow"]
        addfixedstream!(fs, name, flows; ismoleflow=ismoleflow)
    end

    # Create reactions
    reactionnames = config["reactions"]["names"]
    for name in reactionnames
        sname = Symbol(name)
        reactants = config["reactions"][name]["reactants"]
        products = config["reactions"][name]["products"]
        reactantcoeffs = config["reactions"][name]["reactantcoefficients"]
        productcoeffs = config["reactions"][name]["productcoefficients"]
        limiting = config["reactions"][name]["limitingreactant"]

        @eval $(sname)(frac) = Reaction(esc($fs), $reactants, $products, $reactantcoeffs, $productcoeffs, $limiting, frac)
    end

    unitopnames = config["unitops"]["names"]
    for name in unitopnames
        inlets = config["unitops"][name]["inlets"]
        outlets = config["unitops"][name]["outlets"]

        if haskey(config["unitops"][name], "calc")
            calc = Symbol(config["unitops"][name]["calc"])
        else
            calc = nothing
        end

        if haskey(config["unitops"][name], "params")
            params = config["unitops"][name]["params"]
        else
            params = nothing
        end
            
        if isnothing(calc)
            fs.unitops[name] = UnitOp(name, fs.streams, inlets, outlets)
            FlowsheetTools.addtorun!(fs, name)
        else
            @eval $fs.unitops[$name] = UnitOp($name, $fs.streams, $inlets, $outlets, $calc, $params)
            FlowsheetTools.addtorun!(fs, name)
        end
    end

    # Execute all unitops before adding the boundaries, so streams are populated first
    fs()

    boundarynames = config["boundaries"]["names"]
    for name in boundarynames
        unitopnames = config["boundaries"][name]["unitops"]
        fs.boundaries[name] = BalanceBoundary(name, fs.unitops, unitopnames)
    end

    return fs
end
