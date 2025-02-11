using MassBalanceReconciler
using JSON3
using DynamicQuantities

#==============================================================================================
High-level preamble (species, conecntrations etc)
==============================================================================================#
default_massflow = 20.0


clapeyron_pairs = [
    :C1   => "methane",
    :C2   => "ethane",
    :C3   => "propane",
    :iC4  => "isobutane",
    :nC4  => "butane" ,
    :iC5  => "isopentane",
    :nC5  => "pentane",
    :nC6  => "hexane",
    :nC7  => "heptane",
    :nC8  => "octane",
    :nC9  => "nonane",
    :nC10 => "decane",
    :N2   => "nitrogen",
    :CO2  => "carbon dioxide",
    :O2   => "oxygen",
    :H2O  => "water",
    :CO   => "carbon monoxide",
    :H2S  => "hydrogen sulfide",
    :H2   => "hydrogen",
    :He   => "helium",
    :Ar   => "argon",
    :N2O  => "nitrous oxide"
]

LABELS = Tuple( [p[1] for p in clapeyron_pairs] )

default_percent = Species{LABELS}([
    0.4010,
    0.0710,
    0.1030,
    0.0060,
    0.0200,
    0.0000,
    0.0000,
    0.0730,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.3210,
    84.0050,
    0.0000,
    0.0000,
    0.0000,
    15.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000
])

analyzer_labels = Species{LABELS}(split(
"""
Methane
Ethane
Propane
I Butane
N Butane
I Pentane
N Pentane
Hexane
Heptane
Octane
Nonane
Decane
Nitrogen
CO2
Oxygen
H2O
CO
H2S
Hydrogen
Helium
Argon
Nitrous Oxide
""", "\n", keepempty=false
))

default_conc = Species{LABELS}(default_percent ./ sum(default_percent))


#==============================================================================================
Process nodes for a transport system
==============================================================================================#
nodeinfo = NodeInfo[]

push!(nodeinfo, NodeInfo(
    id = :capture,
    stdev   = Dict(LABELS .=> 400.0*900),
    inlets  = [:emissions_source],
    outlets = [:capture_loss, :to_compression]
))

push!(nodeinfo, NodeInfo(
    id = :compression,
    stdev   = Dict(LABELS .=> 400.0*900),
    inlets  = [:to_compression],
    outlets = [:compression_loss, :to_transport] 
))

push!(nodeinfo, NodeInfo(
    id = :transport,
    stdev = Dict(LABELS .=> 400.0*900),
    inlets = [:to_transport],
    outlets = [:transport_loss, :to_injection]
))

push!(nodeinfo, NodeInfo(
    id = :injection,
    stdev = Dict(LABELS .=> 400.0*900),
    inlets  = [:to_injection],
    outlets = [:injection_loss, :to_storage]
))

#==============================================================================================
Process streams for a capture system
==============================================================================================#
streaminfo = StreamInfo[]

push!(streaminfo, StreamInfo(
    id = :emissions_source,
    massflow  = default_massflow,
    molefracs = Dict(default_conc),
    phase = :gas
))

push!(streaminfo, StreamInfo(
    id = :capture_loss,
    massflow = default_massflow*0.001,
    molefracs = :emissions_source,
    phase = :gas
))

push!(streaminfo, StreamInfo(
    id = :to_compression,
    massflow = default_massflow,
    molefracs = :emissions_source,
    phase = :gas
))

push!(streaminfo, StreamInfo(
    id = :compression_loss,
    massflow = default_massflow*0.001,
    molefracs = :to_compression,
    phase = :gas
))

push!(streaminfo, StreamInfo(
    id = :to_transport,
    massflow = default_massflow,
    molefracs = :to_compression,
    phase = :gas
))

push!(streaminfo, StreamInfo(
    id = :transport_loss,
    massflow = default_massflow*0.001,
    molefracs = :to_transport,
    phase = :gas
))

push!(streaminfo, StreamInfo(
    id = :to_injection,
    massflow = default_massflow,
    molefracs = :to_transport,
    phase = :gas
))

push!(streaminfo, StreamInfo(
    id = :injection_loss,
    massflow = default_massflow*0.001,
    molefracs = :to_injection,
    phase = :gas
))

push!(streaminfo, StreamInfo(
    id = :to_storage,
    massflow = default_massflow,
    molefracs = :to_injection,
    phase = :gas
))


#==============================================================================================
Measurements for the capture system
==============================================================================================#
taginfo  = TagInfo[]

push!(taginfo, TagInfo(
    tag   = "emissions_volume_flow",
    units = us"m^3/s",
    stdev = 0.01
))

push!(taginfo, TagInfo(
    tag   = "compressor_mass_flow",
    units = us"kg/s",
    stdev = 0.01
))

push!(taginfo, TagInfo(
    tag   = "transport_mass_flow",
    units = us"kg/s",
    stdev = 0.01
))

push!(taginfo, TagInfo(
    tag   = "injection_mass_flow",
    units = us"kg/s",
    stdev = 0.01
))

for lbl in LABELS
    push!(taginfo, TagInfo(
        tag   = "injection_analyzer."*string(lbl),
        units = Quantity(1.0, SymbolicDimensions()),
        stdev = 0.01
    ))
end

measinfo = MeasInfo[]

push!(measinfo, MeasInfo(
    id = :emissions_volume_flow,
    type = VolumeFlowMeas,
    tags  = Dict(:V=>"emissions_volume_flow", :T=>convert(MeasQuantity, (273.15+25)u"K"), :P=>convert(MeasQuantity, (101.3e3)u"Pa")),
    #stdev = Dict(:V=>0.01),
    stream = :emissions_source
))

push!(measinfo, MeasInfo(
    id = :compressor_mass_flow,
    type = MassFlowMeas,
    tags = Dict(:m=>"compressor_mass_flow"),
    #stdev = Dict(:m=>0.01),
    stream = :to_compression
))

push!(measinfo, MeasInfo(
    id = :transport_mass_flow,
    type = MassFlowMeas,
    tags = Dict(:m=>"transport_mass_flow"),
    #stdev = Dict(:m=>0.01),
    stream = :to_transport
))

push!(measinfo, MeasInfo(
    id = :injection_mass_flow,
    type = MassFlowMeas,
    tags = Dict(:m=>"injection_mass_flow"),
    #stdev = Dict(:m=>0.01),
    stream = :to_injection
))

push!(measinfo, MeasInfo(
    id = :injection_analyzer,
    type  = MoleAnalyzer,
    tags  = Dict(LABELS .=> "injection_analyzer." .* string.(LABELS)),
    #stdev = Dict(LABELS .=> 0.01),
    stream = :to_injection
))

#==============================================================================================
Stream predictions for the system
==============================================================================================#
relationships = StreamRelationship[]

push!(relationships, StreamRelationship(
    id = :injection_loss,
    parent = :to_injection,
    factor = 0.01,
    timeconst = 60.0
))

thermoinfo = ThermoInfo(labels=collect(LABELS), definitions=Dict(clapeyron_pairs))


#==============================================================================================
Assemble entire system
==============================================================================================#
plantinfo = PlantInfo(
    interval = 60.0*15,
    thermo = thermoinfo,
    streams = streaminfo,
    nodes = nodeinfo,
    measurements = measinfo,
    relationships = relationships
)

#==============================================================================================
Output system to JSON
==============================================================================================#
open(joinpath(@__DIR__,"co2_capture_plant.json"), "w") do fh
    JSON3.pretty(fh, plantinfo)
end

capture_dict = open(joinpath(@__DIR__,"co2_capture_plant.json")) do fh
    JSON3.read(fh, numbertype=Float64)
end

reconstructed_plant = PlantInfo(capture_dict)
plantstate = PlantState(reconstructed_plant)
