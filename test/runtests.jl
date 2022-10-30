using FlowsheetTools
using Test

@testset "Streams" begin
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

    @test sysstreams["Test"].moleflows ≈ sysstreams["Product"].moleflows

    copystream!(sysstreams, "Product", "mystream")
    copystream!(sysstreams, "Product", "mystream2"; factor=2.0) # double the flow!

    @test sysstreams["mystream2"].totalmassflow ≈ 2.0 * sysstreams["mystream"].totalmassflow
    @test sysstreams["Product"] !== sysstreams["mystream"]
    @test sysstreams["Product"].atomflows == sysstreams["mystream"].atomflows

    sysstreams["Prod2"] = 2.0*sysstreams["Product"]
    @test sysstreams["Prod2"].totalmassflow ≈ 2.0 * sysstreams["Product"].totalmassflow
end

@testset "UnitOps and Boundaries" begin
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

    sysstreams = StreamList()
    sysunitops = UnitOpList()

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
    
    @boundary begin
        unitops --> ["Reactor", "Membrane"]
    end b sysunitops
    
    # Check the closures
    @test b.atomclosures["C"] ≈ 1.0
    @test b.atomclosures["H"] ≈ 1.0
    @test b.closure ≈ 1.0
    @test b.total_in.totalmassflow ≈ 32.0846
    @test b.total_out.totalmassflow ≈ 32.0846

    @test conversion(b, "Ethane") ≈ 0.0
    @test conversion(b, "Ethylene") ≈ 0.9
    @test selectivity(b, "Ethylene", "Ethane") ≈ 1.0

    copystream!(sysstreams, "Feed", "Feed2"; factor=0.95)
    copystream!(sysstreams, "Product", "Prod2"; factor=1.01)

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
    
    corrections = calccorrections(b2)
    c = round(corrections["Feed2"], sigdigits=6)
    @test c ≈ 1.06043
    c = round(corrections["C2"], sigdigits=6)
    @test c ≈ 1.00731
    c = round(corrections["H2"], sigdigits=6)
    @test c ≈ 1.00791


    b2 = closemb(b2)
    c = round(b2.atomclosures["C"], sigdigits=6)
    @test c ≈ 0.999906
    c = round(b2.atomclosures["H"], sigdigits=6)
    @test c ≈ 1.00007
end
