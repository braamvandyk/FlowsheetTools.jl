using FlowsheetTools
using Test

@testset "Streams" begin
    fs = Flowsheet()
    count = readcomponentlist!(fs, "components", ["Ethylene", "Ethane", "Hydrogen"])
    @test count == 3
    @test fs.comps["Ethylene"].Mr ≈ 28.053
    

    @stream mass begin
        "Ethylene" --> 2.8053
        "Ethane" --> 27.06192
        "Hydrogen" --> 2.21738
    end "Test" fs

    @stream mole begin 
        "Ethylene" --> 0.1
        "Ethane" --> 0.9
        "Hydrogen" --> 1.1
    end "Product" fs

    check = fs.streams["Test"].moleflows .≈ fs.streams["Product"].moleflows
    @test all(values(check))

    copystream!(fs, "Product", "mystream")
    copystream!(fs, "Product", "mystream2"; factor=2.0) # double the flow!

    @test all(values(fs.streams["mystream2"].totalmassflow .≈ 2.0 .* fs.streams["mystream"].totalmassflow))
    @test fs.streams["Product"] !== fs.streams["mystream"]
    @test fs.streams["Product"].atomflows == fs.streams["mystream"].atomflows

    fs.streams["Prod2"] = 2.0*fs.streams["Product"]
    @test all(values(fs.streams["Prod2"].totalmassflow .≈ 2.0 .* fs.streams["Product"].totalmassflow))
end

@testset "UnitOps and Boundaries" begin
    fs = Flowsheet()

    @comp begin
        H --> 2
    end "Hydrogen" fs
    
    @comp begin
        C --> 2
        H --> 6
    end "Ethane" fs
    
    @comp begin
        C --> 2
        H --> 4
    end "Ethylene" fs

    @stream mole begin
        "Ethylene" --> 1.0
        "Hydrogen" --> 2.0
    end "Feed" fs
    
    @stream mole begin
        "Ethylene" --> 0.1
        "Ethane" --> 0.9
        "Hydrogen" --> 1.1
    end "Product" fs
    
    
    @stream mole begin
        "Hydrogen" --> 1.1
    end "H2" fs
    
    @stream mole begin
        "Ethylene" --> 0.1
        "Ethane" --> 0.9
    end "C2" fs
    
    @stream mole begin
        "Hydrogen" --> 0.0
    end "Mixed" fs
    
    @unitop begin
        inlets --> ["Feed"]
        outlets --> ["Product"]
    end "Reactor" fs
    
    @unitop begin
        inlets --> ["Product"]
        outlets --> ["C2", "H2"]
    end "Membrane" fs
    
    @unitop begin
        inlets --> ["H2", "C2"]
        outlets --> ["Mixed"]
        calc --> mixer!
    end "Mixer" fs
    fs.unitops["Mixer"]()
    
    @boundary begin
        unitops --> ["Reactor", "Membrane"]
    end "B" fs
    
    # Check the closures
    @test values(fs.boundaries["B"].atomclosures)[1]["C"] ≈ 1.0
    @test values(fs.boundaries["B"].atomclosures)[1]["H"] ≈ 1.0
    @test values(fs.boundaries["B"].closure)[1] ≈ 1.0
    @test values(fs.boundaries["B"].total_in.totalmassflow)[1] ≈ 32.0846
    @test values(fs.boundaries["B"].total_out.totalmassflow)[1] ≈ 32.0846

    @test conversion(fs.boundaries["B"], "Ethane") ≈ 0.0
    @test conversion(fs.boundaries["B"], "Ethylene") ≈ 0.9
    @test molar_selectivity(fs.boundaries["B"], "Ethylene", "Ethane") ≈ 1.0
end

@testset "Corrections and Closure" begin
    fs = Flowsheet()

    count = readcomponentlist!(fs, "components", ["Ethylene", "Ethane", "Hydrogen", "Nitrogen", "Argon"])
    @test count ==5

    readstreamhistory!(fs, "C2", joinpath("streamhistories", "C2.csv"); ismoleflow=true)
    readstreamhistory!(fs, "H2", joinpath("streamhistories", "Hydrogen.csv"); ismoleflow=true)
    readstreamhistory!(fs, "Product", joinpath("streamhistories", "Product.csv"); ismoleflow=true)
    addemptystream!(fs, "Mixed");
    addemptystream!(fs, "Product1");
    addemptystream!(fs, "Product1a");
    addemptystream!(fs, "Product1b");
    addemptystream!(fs, "Product2");
    addemptystream!(fs, "Product3");
    addemptystream!(fs, "Product4");

    addfixedstream!(fs, "Dummy", [10.0, 0.0, 0.0, 0.0, 0.1])

    @unitop begin
        inlets --> ["H2", "C2"]
        outlets --> ["Mixed"]
        calc --> mixer!
    end "Mixer" fs
    fs.unitops["Mixer"]()

    @test fs.streams["Mixed"] ≈ fs.streams["H2"] + fs.streams["C2"]
    
    @unitop begin
        inlets --> ["Mixed"]
        outlets --> ["Product"]
    end "Reactor" fs

    @unitop begin
        inlets --> ["Product"]
        outlets --> ["Product1", "Product2"]
        calc --> flowsplitter!
        params --> [0.5]
    end "ProductSplitter" fs
    fs.unitops["ProductSplitter"]()

    @test fs.streams["Product"] ≈ fs.streams["Product1"] + fs.streams["Product2"]

    @unitop begin
        inlets --> ["Product1"]
        outlets --> ["Product1a", "Product1b"]
        calc --> componentplitter!
        params --> Dict([
            "Hydrogen" => Dict(["Product1a" => 0.5]),
            "Ethane" => Dict(["Product1b" => 0.3])
        ])
    end "ComponentSplitter" fs
    fs.unitops["ComponentSplitter"]()

    @test fs.streams["Product1"] ≈ fs.streams["Product1a"] + fs.streams["Product1b"]

    @unitop begin
        inlets --> ["Product1a", "Product1b", "Product2"]
        outlets --> ["Product3"]
        calc --> mixer!
    end "Mixer2" fs
    fs.unitops["Mixer2"]()

    @test fs.streams["Product3"] ≈ fs.streams["Product1a"] + fs.streams["Product1b"] + fs.streams["Product2"]

    @test fs.streams["Product"] ≈ fs.streams["Product3"]

    EthaneDehydrogenation(frac) = Reaction(fs, ["Ethane"], ["Ethylene", "Hydrogen"], [1.0], [1.0, 1.0], "Ethane", frac)

    @unitop begin
        inlets --> ["Product3"]
        outlets --> ["Product4"]
        calc --> stoichiometric_reactor!
        params --> [EthaneDehydrogenation(0.5)]
    end "Reactor2" fs
    fs.unitops["Reactor2"]()

    a = values(fs.streams["Product4"].moleflows[Symbol("Ethylene")] .≈ fs.streams["Product4"].moleflows[Symbol("Ethane")])
    @test findall(<(1), a) == [9, 13, 14, 18]

    fs.streams["Mixed"] = 1.1*fs.streams["Mixed"]



    @boundary begin
        unitops --> ["Mixer"]
    end "B1" fs

    @boundary begin
        unitops --> ["Reactor", "ProductSplitter"]
    end "B2" fs

    corrections = calccorrections(fs; λ = 0.0, anchor = "H2")
    @test corrections["C2"] ≈ 1.0
    @test corrections["Product1"] ≈ 1.0
    @test corrections["Mixed"] ≈ 1.0/1.1
    @test corrections["Product2"] ≈ 1.0

    closemb!(fs, corrections)

    sum(values(fs.boundaries["B1"].closure)) ≈ 27
    sum(values(fs.boundaries["B2"].closure)) ≈ 27
end