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
        tags  = ["AI-101-CO2","AI-102-CH4","AI-103-NO2","AI-104-Other"],
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

thermomodel = ThermoModel{ANALYZER_SPECIES}(clapmap)
thermostate = ThermoState{ANALYZER_SPECIES, Float64}(model=thermomodel, T=273+30, P=101.3, n=molepercents./sum(molepercents), phase=:gas)

thermodict = Dict(:s1=>thermostate, :s2=>thermostate, :s3=>thermostate)

#=
#Test conversion/aggregation
mixture = Species{ANALYZER_SPECIES}(rand(length(ANALYZER_SPECIES)))
ghgs = GhgSpecies(mixture)
=#

plantstate = PlantState(plantinfo, thermodict)