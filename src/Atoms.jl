"""
Module Atoms holds the names, symbols and atomic masses of the elements of the Periodic Table.
No exports - use Atoms.atomsymbols and Atoms.atomweights
"""
module Atoms

using DelimitedFiles

_rawatoms = readdlm(joinpath("atoms", "atoms.csv"), ',')
_atomnums = Int.(_rawatoms[:, 1])
_atomsymbols = string.(_rawatoms[:, 2])
_atomweights = Float64.(_rawatoms[:, 4])
atomsymbols = Dict(_atomnums .=> _atomsymbols)
atomweights = Dict(_atomsymbols .=> _atomweights)

end
