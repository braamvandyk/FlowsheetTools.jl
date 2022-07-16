#----Component---------------

struct Component
    name
    atoms
    counts
    Mr::Float64
end
StructTypes.StructType(::Type{Component}) = StructTypes.Struct()


"""
    Component(name, atoms, counts)

Constructor for molecular species that defines the atomic make-up.
"""
function Component(name, atoms, counts)
    Mr = 0.0
    for (i, atom) in enumerate(atoms)
        Mr += counts[i]*atomweights[atom]
    end
    Component(name, atoms, counts, Mr)
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


# Parsing an expression to extract atoms and counts for use in macro @comp
function parse_comp(ex, name)
    atoms = String[]
    counts = Int[]
    
    for line in ex.args
        match_atom = @capture(line, atom_sym_ --> count_)
        if match_atom
            atom = String(atom_sym)
            if !(atom in values(atomsymbols))
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
    return :(comp = Component($name, $atoms, $counts))
end


"""
    @comp begin
        C --> 2
        H --> 6
    end "Ethane"

Defines a Component() with the specified name and atomic composition.
"""
macro comp(ex::Expr, name::String)      
    return parse_comp(ex, name)
end


"""
    function writecomp(filename, comp)

Write a Component struct to a JSON file.
"""
function writecomponent(comp)
    filename = lowercase(joinpath("components", comp.name * ".json"))
    open(filename, "w") do io
        JSON3.pretty(io, comp)
    end
end


"""
    function readcomp(filename)

Read a Component struct from a JSON file.
"""
function readcomponent(filename)
    open(filename, "r") do io
        JSON3.read(io, Component)
    end
end

function readcomponentlist(filenames)
    count = 0
    res = Component[]
    available = readdir("components/")
    
    for fn in filenames
        fname = lowercase(fn * ".json")
        if fname in available
            comp = readcomponent(joinpath("components", fname))
            push!(res, comp)
            count += 1
        end
    end

    return res, count
end