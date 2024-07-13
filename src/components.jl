#----------------------------------------------------------------------------
# 
#----Definitions-------------------------------------------------------------
# 
#----------------------------------------------------------------------------

struct Component
    name::String
    atoms::Vector{String}
    counts::Vector{Int}
    Mr::Float64
end

struct ComponentList
    list::OrderedDict{String, Component}
end


#----------------------------------------------------------------------------
# 
#----Constructors------------------------------------------------------------
# 
#----------------------------------------------------------------------------


"""

    Component(name, atoms, counts)

Constructor for molecular species that defines the atomic make-up.
It is recommended to rather use the @comp macro to create components interactively.
The atoms are identified by their IUPAC symbols, as a `String`. Counts are integers.

"""
function Component(name, atoms, counts)
    Mr = 0.0
    for (i, atom) in enumerate(atoms)
        Mr += counts[i]*Atoms.atomweights[atom]
    end
    return Component(name, atoms, counts, Mr)
end


"""

    ComponentList()

Constructor for an empty component list. Components are added when created via the @comp macro.

"""
function ComponentList()
    l = OrderedDict{String, Component}()
    return ComponentList(l)
end


#----------------------------------------------------------------------------
#
#----Base overloads----------------------------------------------------------
#
#----------------------------------------------------------------------------


function Base.setindex!(A::ComponentList, X::Component, idx::String)
    A.list[idx] = X

    return nothing
end

function Base.getindex(A::ComponentList, idx::String)
    return A.list[idx]
end

function Base.length(A::ComponentList)
    return length(A.list)
end

function Base.iterate(A::ComponentList)
    return iterate(A.list)
end

function Base.iterate(A::ComponentList, state)
    return iterate(A.list, state)
end


#----------------------------------------------------------------------------
#
#----Pretty Printing---------------------------------------------------------
#
#----------------------------------------------------------------------------


function Base.show(io::IO, c::Component)
    println(io, "Component: $(c.name)\n")
    println(io, "  Atom\t\tCount")
    println(io, "-"^21)
    for (i, atom) in enumerate(c.atoms)
        println(io, "   ", rpad(atom,2), "\t\t", lpad(c.counts[i], 4))
    end
end


function  Base.show(io::IO, cl::ComponentList)
    println(io, "Component list:")
    for comp in cl.list
        println(io, "  ", comp.first)
    end
end


#----------------------------------------------------------------------------
#
#----Macros------------------------------------------------------------------
#
#----------------------------------------------------------------------------


"""

    @comp begin
        C --> 2
        H --> 6
    end "Ethane" fs

Defines a Component with the specified name and atomic composition and add it to fs.comps, a ComponentList, inside fs::Flowsheet.

"""
macro comp(ex::Expr, name::String, fs::Symbol)      
    local atoms = String[]
    local counts = Int[]

    for line in ex.args
        match_atom = @capture(line, atom_sym_ --> count_)
        if match_atom
            atom = String(atom_sym)
            if !(atom in values(Atoms.atomsymbols))
                throw(ArgumentError("Invalid symbol for atom specified."))
            end
            if atom in atoms
                i = findfirst(x -> x == atom, atoms)
                counts[i] += count
            else
                push!(atoms, atom)
                push!(counts, count)
            end
        end
    end
    
    return :($(esc(fs)).comps[$name] = Component($name, $atoms, $counts))
end


#----------------------------------------------------------------------------
#
#----Utilities---------------------------------------------------------------
#
#----------------------------------------------------------------------------


"""

    writecomponent(filename, comp)

Write a Component struct to a file.

"""
function writecomponent(filename, comp)
    str = "$(comp.name)\n"
    str = str*string(length(comp.atoms))*"\n"
    for atom in comp.atoms
        str = str*atom*"\n"
    end
    for count in comp.counts
        str = str*string(count)*"\n"
    end
    str = str*string(comp.Mr)
    
    filename = lowercase(filename)   
    open(filename, "w") do io
        write(io, str)
    end
end


"""

*Internal use only!* Use readcomponentlist!() instead.

    function readcomp(filename)

Read a Component struct from a file. This file is created using `writecomponent()` for a Component object, created via @comp.

"""
function readcomponent(filename)
    f = open(filename, "r")
    lines = readlines(f)
    close(f)

    name = lines[1]
    count = parse(Int, lines[2])

    atoms = Vector{String}(undef, count)
    for i in 1:count
        atoms[i] = lines[i+2]
    end

    counts = Vector{Int}(undef, count)
    for i in 1:count
        counts[i] = parse(Int, lines[i+2+count])
    end

    Mr = parse(Float64, lines[end])

    return Component(name, atoms, counts, Mr)
end


"""

    function readcomponentlist!(fs, foldername, filenames)

Read a list of components into a ComponentList, fs.comps.
The folder is specified and `filenames` contains a list of filenames without extentions.
The component files are text files, created with `writecomponent()` for a Component object, created via @comp.

"""
function readcomponentlist!(fs, foldername, filenames)
    @argcheck fs isa Flowsheet "fs must be a Flowsheet"
    count = 0
    # Get the list of files available in the folder
    available = readdir(foldername)
    
    for fn in filenames
        fname = fn * ".comp"
        if lowercase(fname) in available
            fs.comps[fn] = readcomponent(joinpath(foldername, fname))
            count += 1
        else
            @warn "Component file $(joinpath(foldername, fname)) not found."
        end
    end

    return count
end


"""

    function deletecomponent!(fs, name)

Delete a component from fs::Flowsheet's ComponentList. Identified via name::String.
This will result in all streams and mass balance boundaries also being deleted, as 
these all refer to all components.

"""
function deletecomponent!(fs, name)
    @argcheck fs isa Flowsheet "fs must be a Flowsheet"
    if name in keys(fs.comps.list)
        delete!(fs.comps.list, name)
    end

    # Now also clear all stream, unitops and boundaries
    deletestreams!(fs)
end


"""

    function deletecomponents!(fs)

Delete all components from the Flowsheet's ComponentList.
This will result in all streams and mass balance boundaries also being deleted, as 
these all refer to all components.

"""
function deletecomponents!(fs)
    @argcheck fs isa Flowsheet "fs must be a Flowsheet"
    for name in keys(fs.comps.list)
        delete!(fs.comps.list, name)
    end

    # Now also clear all stream, unitops and boundaries
    deletestreams!(fs)

end