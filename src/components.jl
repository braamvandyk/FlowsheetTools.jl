#----------------------------------------------------------------------------
# 
#----Definitions-------------------------------------------------------------
# 
#----------------------------------------------------------------------------


struct Component
    name
    atoms
    counts
    Mr::Float64
end
StructTypes.StructType(::Type{Component}) = StructTypes.Struct()


struct ComponentList
    list
end


#----------------------------------------------------------------------------
# 
#----Constructors------------------------------------------------------------
# 
#----------------------------------------------------------------------------


"""
    Component(name, atoms, counts)

Constructor for molecular species that defines the atomic make-up.
It is recommended to rather use the @comp macro to create components.
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

Constructor for an empty component list. Components are added when created via the @comp macro
"""
function ComponentList()
    l = Dict{String, Component}()
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


# Pretty printing for component objects
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
    end "Ethane" syscomps

Defines a Component with the specified name and atomic composition and add it to syscomps::ComponenList
"""
macro comp(ex::Expr, name::String, complist::Symbol)      
    local atoms = String[]
    local counts = Int[]

    for line in ex.args
        match_atom = @capture(line, atom_sym_ --> count_)
        if match_atom
            atom = String(atom_sym)
            if !(atom in values(Atoms.atomsymbols))
                error("Invalid symbol for atom specified.")
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
    
    return :($(esc(complist))[$name] = Component($name, $atoms, $counts))
end


#----------------------------------------------------------------------------
#
#----Utilities---------------------------------------------------------------
#
#----------------------------------------------------------------------------


"""
    function writecomp(filename, comp)

Write a Component struct to a JSON file.
"""
function writecomponent(filename, comp)
    filename = lowercase(filename)
    open(filename, "w") do io
        JSON3.pretty(io, comp)
    end
end


"""
    function readcomp(filename::String)

Read a Component struct from a JSON file.
"""
function readcomponent(filename::String)
    open(filename, "r") do io
        JSON3.read(io, Component)
    end
end


"""
    function readcomponentlist(complist::ComponentList, folder::String, filenames::Vector{String})

Read a list of components into a ComponentList
The folder is specified and `filenames` contains a list of filenames without extentions.
The component files are .json files created with `writecomp()`
"""
function readcomponentlist!(complist::ComponentList, folder::String, filenames::Vector{String})
    count = 0
    # Get the list of files available in the folder
    available = readdir(folder)
    
    for fn in filenames
        fname = fn * ".json"
        if lowercase(fname) in available
            complist[fn] = readcomponent(joinpath(folder, fname))
            count += 1
        end
    end

    return count
end