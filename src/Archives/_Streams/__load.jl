
include(joinpath(@__DIR__,"Chem.jl"))
include(joinpath(@__DIR__,"Streams.jl"))
include(joinpath(@__DIR__,"StreamThermo.jl"))
include(joinpath(@__DIR__,"StreamDispatch.jl"))

CHEM = compile_chemical_data()
build_component_indexer(Tuple(CHEM.Names))