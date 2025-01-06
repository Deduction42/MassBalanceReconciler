using MassBalanceReconciler
using JSON3

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
nodes = NodeInfo[]

push!(nodes, NodeInfo(
    id = :capture,
    stdev   = Dict(LABELS .=> 30.0),
    inlets  = [:emissions_source],
    outlets = [:capture_loss, :to_compression]
))

push!(nodes, NodeInfo(
    id = :compression,
    stdev   = Dict(LABELS .=> 30.0),
    inlets  = [:to_compression],
    outlets = [:compression_loss, :to_transport] 
))

push!(nodes, NodeInfo(
    id = :transport,
    stdev = Dict(LABELS .=> 30.0),
    inlets = [:to_transport],
    outlets = [:transport_loss, :to_injection]
))

push!(nodes, NodeInfo(
    id = :injection,
    stdev = Dict(LABELS .=> 30.0),
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
