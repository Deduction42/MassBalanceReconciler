include("_Species.jl")

using JSON3
using FlexUnits, .UnitRegistry
using JSON3.StructTypes
using Accessors

import FlexUnits.dimension, FlexUnits.ubase
import Base.RefValue
import Base.Fix1
import Base.Fix2

#Settings for default quantities
abstract type AbstractInfo end
const MeasUnits = AffineUnits{FlexUnits.DEFAULT_DIMENSONS}
const MeasQuantity = Quantity{Float64, FlexUnits.DEFAULT_DIMENSONS}


const NULL_SYMBOL = :_
const MOL_ε = 1e-9

#=======================================================================================
# Structure that contains thermodynamic information to build thermo models
=======================================================================================#
@kwdef struct ThermoInfo
    labels  :: Vector{Symbol}
    definitions :: Dict{Symbol, Union{String, Dict{String,Float64}}}
end

function ThermoInfo(labels::AbstractVector, definitions::AbstractDict)
    T = Union{String, Dict{String,Float64}}
    return ThermoInfo(
        labels = symbolize(labels),
        definitions = Dict{Symbol, T}(Symbol(k)=>_thermodef(v) for (k,v) in pairs(definitions))
    )
end

function ThermoInfo(d::AbstractDict{Symbol})
    return ThermoInfo(
        labels = symbolize(d[:labels]),
        definitions = d[:definitions]
    )
end

ThermoInfo(d::AbstractDict{<:AbstractString}) = symbolize(d)


#=============================================================================
Construction info for streams
=============================================================================#
@kwdef struct StreamInfo <: AbstractInfo
    id        :: Symbol
    massflow  :: Float64
    molefracs :: Dict{Symbol, Float64}
    copycomp  :: RefValue{Symbol} = Ref(NULL_SYMBOL)
    phase     :: Symbol = :unknown
end

massflow(info::StreamInfo)  = info.massflow 
molefracs(info::StreamInfo) = info.molefracs

function StreamInfo(d::AbstractDict{Symbol})
    molefracs = d[:molefracs]

    return StreamInfo(
        id = Symbol(d[:id]),
        massflow  = d[:massflow],
        molefracs = (molefracs isa AbstractDict) ? symbolize(Float64, molefracs) : symbolize(molefracs),
        phase = Symbol(get(d, :phase, :unknown))
    )
end

StreamInfo(d::AbstractDict{<:AbstractString}) = StreamInfo(symbolize(d))

@kwdef struct TagInfo <: AbstractInfo
    name  :: String = ""
    val   :: Float64 = NaN
    units :: MeasUnits = u""
end
getname(t::TagInfo)  = t.name 
getvalue(t::TagInfo) = t.val 
getunits(t::TagInfo) = t.units

function TagInfo(d::AbstractDict{Symbol}) 
    #Convert units to SI if there are constant values
    rawtag = d[:tag]
    name = (rawtag isa AbstractString) ? String(rawtag) : ""
    val = (rawtag isa Real) ? Float64(rawtag) : NaN
    units = parse_units(d[:units])
    return ubase( TagInfo(name=name, val=val, units=units) )
end

dimension(tag::TagInfo) = dimension(tag.units)
hasname(tag::TagInfo) = !isempty(tag.name)

function ubase(tag::TagInfo) #Convert the value to SI units
    if hasname(tag)
        return tag 
    else
        basequant = ubase(tag.val*tag.units)
        return TagInfo(name=tag.name, val=ustrip(basequant), units=unit(basequant))
    end
end

#=============================================================================
Construction info for nodes
=============================================================================#

@kwdef struct NodeInfo <: AbstractInfo
    id        :: Symbol
    stdev     :: Float64
    inlets    :: Vector{Symbol}
    outlets   :: Vector{Symbol}
    reactions :: Vector{Dict{Symbol, Float64}} = Dict{Symbol, Float64}[]
    uniform_split :: Bool = false
end

function NodeInfo(d::AbstractDict{Symbol})
    nodeinfo = NodeInfo(
        id = Symbol(d[:id]),
        stdev  = d[:stdev],
        inlets = symbolize(d[:inlets]),
        outlets = symbolize(d[:outlets]),
        reactions = symbolize.(Float64, d[:reactions])
    )

    !isempty(nodeinfo.inlets)  || error("Error parsing node '$(nodeinfo.id)': Inlets cannot be empty")
    !isempty(nodeinfo.outlets) || error("Error parsing node '$(nodeinfo.id)': Outlets cannot be empty")

    return nodeinfo
end

NodeInfo(d::AbstractDict{<:AbstractString}) = NodeInfo(symbolize(d))

function add_reaction!(nodeinfo::NodeInfo, stoich::Species)
    push!(nodeinfo.reactions, ReactionRef{L,N}(0, stoich))
end

function composition_is_conserved(n::NodeInfo)
    return isempty(n.reactions) && (length(n.inlets)<=1 || n.uniform_split)
end

#=============================================================================
Fill copycomps and connetion in StreamInfos using NodeInfos
=============================================================================#
function fill_copycomps!(streams::AbstractVector{StreamInfo}, nodes::AbstractVector{NodeInfo})
    stream_dict = Dict( stream.id=>stream for stream in streams )
    id_massflow(id::Symbol) = massflow(stream_dict[id])

    for node in Iterators.filter(composition_is_conserved, nodes)
        sort!(node.outlets, rev=true, by=id_massflow) #Ensure the outlets are sorted by greatest massflow first

        if length(node.inlets) == 1 #Pass inlet compositions to the outlet
            copy_id = node.inlets[1]

        elseif length(node.inlets) > 1 #Set the first (largest) outlet as the reference, copy its composition to other outlets
            copy_id = node.outlets[1]

        else
            error("Number of inlets must be greater than zero")
        end

        for outlet in node.outlets
            if outlet == copy_id #Don't create circular references
                stream_dict[outlet].copycomp[] = NULL_SYMBOL
            else
                stream_dict[outlet].copycomp[] = copy_id
            end
        end
    end

    return streams
end


#==========================================================================================================
Construction info for measurements
==========================================================================================================#
@kwdef struct MeasInfo <: AbstractInfo
    id     :: Symbol
    type   :: UnionAll
    tags   :: Dict{Symbol, TagInfo}
    stdev  :: Quantity
    stream :: Symbol = NULL_SYMBOL
    node   :: Symbol = NULL_SYMBOL
end

function MeasInfo(d::AbstractDict{Symbol}; stream=NULL_SYMBOL, node=NULL_SYMBOL) 
    MeasType = eval(Meta.parse(d[:type]))

    function buildtags(tagd)
        return Dict{Symbol, TagInfo}(Symbol(k)=>TagInfo(v) for (k, v) in pairs(tagd))
    end

    return MeasInfo(
        id     = Symbol(d[:id]),
        type   = MeasType,
        tags   = buildtags(d[:tags]),
        stdev  = parse_stdev(d[:stdev]),
        stream = Symbol(get(d, :stream, stream)),
        node   = Symbol(get(d, :node, node))
    )
end

MeasInfo(d::AbstractDict{<:AbstractString}) = MeasInfo(symbolize(d))


#=============================================================================
Construction info for entire system
=============================================================================#
@kwdef struct PlantInfo <: AbstractInfo
    interval :: Float64
    thermo   :: ThermoInfo
    streams  :: Vector{StreamInfo} = StreamInfo[]
    nodes    :: Vector{NodeInfo}   = NodeInfo[]
    measurements  :: Vector{MeasInfo}  = MeasInfo[]
end

function PlantInfo(d::AbstractDict{<:Symbol})
    return PlantInfo(
        interval = d[:interval],
        thermo  = ThermoInfo(d[:thermo]),
        streams = StreamInfo.(d[:streams]),
        nodes   = NodeInfo.(d[:nodes]),
        measurements = MeasInfo.(d[:measurements])
    )
end

function fill_copycomps!(plantinfo::PlantInfo) 
    fill_copycomps!(plantinfo.streams, plantinfo.nodes)
    return plantinfo
end

#==========================================================================================================
Utility functions
==========================================================================================================#
StructTypes.StructType(::Type{MeasQuantity}) = StructTypes.StringType()
StructTypes.StructType(::Type{MeasUnits}) = StructTypes.StringType()

MeasUnits(x::String) = parse_units(x)
MeasQuantity(x::String) = parse_quantity(x)


function parse_stdev(x::AbstractString)
    stdev = parse_quantity(x)
    if !iszero(FlexUnits.uoffset(unit(stdev)))
        return FlexUnits.remove_offset(stdev)
    else
        return stdev 
    end
end

function parse_quantity(x::AbstractString)
    return qparse(x)
end

function parse_units(x::AbstractString)
    return uparse(x)
end


function tryparse_quantity(x::AbstractString)
    try 
        return parse_quantity(x::AbstractString)
    catch
        return String(x)
    end
end

symbolize(x::AbstractString) = Symbol(x)
symbolize(x::AbstractVector{<:AbstractString}) = Symbol.(x)
symbolize(x::AbstractDict{<:AbstractString, T}) where T = Dict(Symbol(k)=>v for (k,v) in pairs(x))
symbolize(x::AbstractDict{<:AbstractString, <:Real}) = Dict(Symbol(k)=>Float64(v) for (k,v) in pairs(x))

symbolize(x::Symbol) = x
symbolize(x::AbstractVector{Symbol}) = Vector(x) 
symbolize(x::AbstractDict{Symbol, T}) where T = Dict(Symbol(k)=>v for (k,v) in pairs(x))
symbolize(x::AbstractDict{Symbol, <:Real}) = Dict(Symbol(k)=>Float64(v) for (k,v) in pairs(x))

symbolize(::Type{T}, x::AbstractDict) where T = Dict{Symbol, T}(Symbol(k)=>convert(T, v) for (k,v) in pairs(x))
#symbolize(::Type{T}, x::AbstractDict) where T <: AbstractQuantity = Dict{Symbol, Union{String, MeasQuantity}}(Symbol(k)=>tryparse_units(v) for (k,v) in pairs(x))
