using MassBalanceReconciler
using JSON3
using Revise


#include("_GhgSpecies.jl")

function label2symbol(name::AbstractString)
    return Symbol(replace(lowercase(name), r"\s+"=>"_"))
end

clappairs = [
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
    :argon => "argon",
    :n2o => "nitrous oxide"
]
clapmap = Dict(clappairs)
GASES = Tuple(p[1] for p in clappairs)

analyzerjson = open(joinpath(@__DIR__,"analyzer.json")) do fh 
    JSON3.read(fh)
end

analyzerdict = Dict(label2symbol.(analyzerjson.components) .=> analyzerjson.mole_percents)
analyzerspec = Species{GASES}(analyzerdict)
moletags = Species{GASES}("AI-101 ".*analyzerjson.components)

streamcomp = begin
    v = analyzerspec[:] .+ 0.1
    Species{GASES}(v./sum(v))
end
compdict = Dict{Symbol, Float64}(GASES .=> streamcomp[:])

nodeinfo = [
    NodeInfo(
        id = :v1,
        stdev = Dict(GASES .=> 10.0),
        inlets = [:s1, :s2],
        outlets = [:s3]
    )
]

streaminfo = [
    StreamInfo(
        id = :s1,
        massflow = 20.0,
        molefracs = compdict,
        phase = :gas
    ),
    StreamInfo(
        id = :s2,
        massflow = 10.0,
        molefracs = compdict,
        phase = :gas
    ),
    StreamInfo(
        id = :s3,
        massflow = 10.0,
        molefracs = compdict,
        phase = :gas
    )
]

measinfo = [
    MeasInfo(
        id = :s1_volume,
        type  = VolumeFlowMeas,
        tags  = Dict(:V =>"FI-101", :T =>273.15+100, :P =>10*101.3e3),
        stdev = Dict(:V => 0.01),
        stream = :s1
    ),
    MeasInfo(
        id = :s2_mass,
        type  = MassFlowMeas,
        tags  = Dict(:m =>"FI-102"),
        stdev = Dict(:m =>0.001),
        stream = :s2
    ),
    MeasInfo(
        id    = :s1_analyzer,
        type  = MoleAnalyzer,
        tags  = Dict(GASES .=> moletags),
        stdev = Dict(GASES .=> 0.01),
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

plantinfo = PlantInfo(
    species=collect(GASES), 
    thermo=clapmap,
    streams = streaminfo,
    nodes = nodeinfo,
    measurements = measinfo,
    relationships = relationships
)



plantstate = PlantState(plantinfo)



tagdict = Dict{String,Float64}()
tagdict["TI-101"] = 273.15+35
tagdict["PI-101"] = 104.3
tagdict["FI-101"] = 30.0
tagdict["FI-102"] = 13.0

fracs = rand(length(GASES))
fracs = fracs./sum(fracs)
for (ii, k) in enumerate("AI-101 ".*analyzerjson.components)
    tagdict[k] = fracs[ii]
end


readvalues!(plantstate, tagdict, 60.0*15)
updatethermo!(plantstate)
predict!(plantstate, 60)


statevec = deepcopy(plantstate.statevec)
negloglik(plantstate.statevec, plantstate)

@time results = reconcile_statevec!(plantstate)
P0 = deepcopy(plantstate.statecov)
@time P1 = reconcile_statecov!(plantstate)
display(P0)
display(P1)