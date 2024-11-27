using MassBalanceReconciler
using Revise
using JSON3

include("_GhgSpecies.jl")

function label2symbol(name::AbstractString)
    return Symbol(replace(lowercase(name), r"\s+"=>"_"))
end

clapmap = Dict{Symbol,String}(
    :methane => "methane",
    :ethane => "ethane",
    :propane => "propane",
    :i_butane => "isobutane",
    :n_butane => "butane" ,
    :i_pentane => "isopentane",
    :n_pentane => "pentane",
    :hexane => "hexane",
    :heptane => "heptane",
    :octane  => "octane",
    :nonane => "nonane",
    :decane => "decane",
    :nitrogen => "nitrogen",
    :co2 => "carbon dioxide",
    :oxygen => "oxygen",
    :h2o => "water",
    :co => "carbon monoxide",
    :h2s => "hydrogen sulfide",
    :hydrogen => "hydrogen",
    :helium => "helium",
    :argon => "argon"
)

GASES = Tuple(collect(keys(clapmap)))

plantinfo = PlantInfo{GHG}()

nodeinfo = [
    NodeInfo{GHG}(
        id = :v1,
        stdev = Species{GHG}([10.0, 10.0, 10.0, 20.0]),
        inlets = [:s1, :s2],
        outlets = [:s3]
    )
]

streaminfo = [
    StreamInfo{GHG}(
        id = :s1,
        massflow = 20.0,
        phase = :gas
    ),
    StreamInfo{GHG}(
        id = :s2,
        massflow = 10.0,
        phase = :gas
    ),
    StreamInfo{GHG}(
        id = :s3,
        massflow = 10.0,
        phase = :gas
    )
]

measinfo = [
    MeasInfo(
        id = :s1_volume,
        type  = VolumeFlowMeas,
        tags  = ["FI-101"],
        stdev = [0.1],
        stream = :s1
    ),
    MeasInfo(
        id = :s2_mass,
        type = MassFlowMeas,
        tags = ["FI-102"],
        stdev = [0.01],
        stream = :s2
    ),
    MeasInfo(
        id    = :s1_analyzer,
        type  = MoleAnalyzer,
        tags  = ["AI-101 CO2","AI-101 CH4","AI-101 NO2","AI-101 Other"],
        stdev = [0.01, 0.01, 0.01, 0.01],
        stream = :s1
    )
]

relationships = [
    StreamRelationship(
        id = :s2,
        parent = :s1,
        factor = 0.5,
        timeconst = 60.0
    )
]

plantinfo = PlantInfo{GHG}(
    streams = streaminfo,
    nodes = nodeinfo,
    measurements = measinfo,
    relationships = relationships
)

jsonobj = open(joinpath(@__DIR__,"analyzer.json")) do fh 
    JSON3.read(fh)
end

const ANALYZER_SPECIES = Tuple(label2symbol.(jsonobj.components))
const GHG_MAP = GhgSpecies{Union{Symbol,Nothing}}(
    CO2   = :co2,
    CH4   = :methane,
    N2O   = nothing,
    other = nothing
)
GhgSpecies(mixture::Species{ANALYZER_SPECIES}) = GhgSpecies(mixture, GHG_MAP)
GhgSpecies{T}(mixture::Species{ANALYZER_SPECIES}) where T = GhgSpecies{T}(mixture, GHG_MAP)

MassBalanceReconciler.Species{GHG}(mixture::Species{ANALYZER_SPECIES}) = GhgSpecies(mixture)
MassBalanceReconciler.Species{GHG,T}(mixture::Species{ANALYZER_SPECIES}) where T = GhgSpecies{T}(mixture)

moledict = Dict(label2symbol.(jsonobj.components) .=> jsonobj.mole_percents)
molepercents = Species{ANALYZER_SPECIES}([moledict[k] for k in ANALYZER_SPECIES])
moletags = Species{ANALYZER_SPECIES}("AI-101 ".*jsonobj.components)

thermomodel = ThermoModel{ANALYZER_SPECIES}(clapmap)
thermostate = ThermoState{ANALYZER_SPECIES, Float64}(model=thermomodel, T=273+30, P=101.3, n=molepercents./sum(molepercents), phase=:gas)
thermotags  = ThermoState{ANALYZER_SPECIES, String}(model=thermomodel, T="TI-101", P="PI-101", n=moletags, phase=:gas)


thermoinfo = ThermoInfo{ANALYZER_SPECIES}(
    tags   = Dict(:s1=>thermotags, :s2=>thermotags, :s3=>thermotags),
    values = Dict(:s1=>thermostate, :s2=>thermostate, :s3=>thermostate)
)




#=
#Test conversion/aggregation
mixture = Species{ANALYZER_SPECIES}(rand(length(ANALYZER_SPECIES)))
ghgs = GhgSpecies(mixture)
=#

plantstate = PlantState(plantinfo, thermoinfo.values)



tagdict = Dict{String,Float64}()
tagdict["TI-101"] = 273.15+35
tagdict["PI-101"] = 104.3
tagdict["FI-101"] = 30.0
tagdict["FI-102"] = 13.0

fracs = rand(length(ANALYZER_SPECIES))
fracs = fracs./sum(fracs)
for (ii, k) in enumerate("AI-101 ".*jsonobj.components)
    tagdict[k] = fracs[ii]
end

readvalues!(thermoinfo, tagdict)
translate!(tagdict, plantstate.measurements, thermoinfo.values)
updatethermo!(plantstate, thermoinfo.values)
readvalues!(plantstate, tagdict)
