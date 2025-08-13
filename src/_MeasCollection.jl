include("_AbstractAnalyzer.jl")
abstract type AbstractMeasCollection{S,T,N} end

basetype(T::DataType) = T.name.wrapper
basetype(T::UnionAll) = basetype(T.body)

#==========================================================================================================
Mole Balance around each node (counts as observation)
==========================================================================================================#
@kwdef mutable struct MoleBalance{S, T, N} <: AbstractMultiMeas{S, T}
    id        :: Symbol
    value     :: Species{S, T, N}
    interval  :: Float64
    inlets    :: Vector{StreamRef{S, N}}
    outlets   :: Vector{StreamRef{S, N}}
    reactions :: Vector{ReactionRef{S, N}}
    stdev     :: Float64
end
MoleBalance{S, T}(x...) where {S,T} = MoleBalance{S, T, length(S)}(x...)
MoleBalance{S, T}(;kw...) where {S,T} = MoleBalance{S, T, length(S)}(;kw...)

function MoleBalance(nodeinfo::NodeInfo, streams::Dict{Symbol, <:StreamRef{S}}, interval::Base.RefValue{Float64}) where S
    nodeid = nodeinfo.id
    N = length(S)

    return MoleBalance{S, Float64}(
        id        = nodeid,
        value     = zero(Species{S, Float64, N}),
        interval  = interval[],
        stdev     = nodeinfo.stdev,
        inlets    = [streams[id] for id in nodeinfo.inlets],
        outlets   = [streams[id] for id in nodeinfo.outlets],
        reactions = nodeinfo.reactions
    )
end

#Mole balances don't get read from measurement dictionaries, they are propagated from state results
readvalues!(m::MoleBalance, d::Dict) = m

function addinds!(inds::BitVector, m::MoleBalance)
    addinds!(inds, m.inlets)
    addinds!(inds, m.outlets)
    addinds!(inds, m.reactions)
    return inds
end

function prediction(x::AbstractVector{T}, m::MoleBalance{S, <:Float64, N}) where {S,T,N}
    RT = promote_type(T,Float64)
    default = zero(SVector{N,RT})

    balance = (
          (isempty(m.inlets)    ? default : sum(Base.Fix1(speciesvec, x), m.inlets))
        - (isempty(m.outlets)   ? default : sum(Base.Fix1(speciesvec, x), m.outlets))
        + (isempty(m.reactions) ? default : sum(Base.Fix1(speciesvec, x), m.reactions))
    )

    return balance.*m.interval
end

function loglik(x::AbstractVector{T}, m::MoleBalance) where T <: Real
    return -0.5*sum(x->abs2(x/getstdev(m)), innovation(x, m))
end

getvariance(m::MoleBalance{S}) where S = @SArray fill(abs2(m.stdev), length(S))


#=============================================================================
Collection of all measurements
=============================================================================#
@kwdef struct MeasCollection{S,T,N} <: AbstractMeasCollection{S,T,N}
    MassFlowMeas    :: Vector{MassFlowMeas{S,T,N}}   = MassFlowMeas{S,T,N}[]
    VolumeDensMeas  :: Vector{VolumeDensMeas{S,T,N}} = VolumeDensMeas{S,T,N}[]
    VolumeFlowMeas  :: Vector{VolumeFlowMeas{S,T,N}} = VolumeFlowMeas{S,T,N}[]
    MoleAnalyzer    :: Vector{MoleAnalyzer{S,T,N}}   = MoleAnalyzer{S,T,N}[]
    MassAnalyzer    :: Vector{MassAnalyzer{S,T,N}}   = MassAnalyzer{S,T,N}[]
    MoleBalance     :: Vector{MoleBalance{S,T,N}}    = MoleBalance{S,T,N}[]
end
MeasCollection{S,T}(;kwargs...) where {S,T} = MeasCollection{S,T,length(S)}(kwargs...)

Base.getindex(m::MC, k::Symbol) where MC<:AbstractMeasCollection  = getproperty(m, k)
Base.getindex(m::MC, k::Integer) where MC<:AbstractMeasCollection = getproperty(m, fieldnames(MC)[k])
Base.getindex(m::MC, k::Colon) where MC<:AbstractMeasCollection  = map(Base.Fix1(getproperty, m), fieldnames(MC))
Base.firstindex(m::MC) where MC<:AbstractMeasCollection = 1
Base.lastindex(m::MC) where MC<:AbstractMeasCollection = length(fieldnames(MC))
Base.getindex(m::MC, ::Type{T}) where {T, MC<:AbstractMeasCollection} = getproperty(m, Symbol(basetype(T)))
Base.getindex(m::MC, k::AbstractVector) where MC<:AbstractMeasCollection = map(Base.Fix1(getindex, m), k)
Base.getindex(m::MC, k::Tuple) where MC<:AbstractMeasCollection = map(Base.Fix1(getindex, m), k)

function Base.getindex(m::MC, k::MeasReference) where MC<:AbstractMeasCollection
    return m[k.type][k]
end

function populate!(c::AbstractMeasCollection, info::MeasInfo, streams::Dict{Symbol, StreamRef})
    type = info.type
    return push!(c[type], type(info, streams))
end

function readvalues!(c::MC, d::Dict, ind::Integer...) where MC<:AbstractMeasCollection
    for fn in fieldnames(MeasCollection)
        readvalues!(c[fn], d, ind...)
    end
    return c
end

function readvalues!(vmeas::AbstractVector{M}, d::Dict, ind::Integer...) where {M <: AbstractMeas}
    for m in vmeas 
        readvalues!(m, d, ind...)
    end
    return vmeas
end

function updatethermo!(vmeas::AbstractVector{M}, statevec::AbstractVector{<:Real}, thermo::ThermoModel) where {M <: AbstractMeas}
    for m in vmeas
        updatethermo!(m, statevec, thermo)
    end
    return vmeas
end

function updatethermo!(c::MC, statevec::AbstractVector{<:Real}, thermo::ThermoModel) where MC <: AbstractMeasCollection
    map(fn -> updatethermo!(c[fn], statevec, thermo), fieldnames(MC))
    return c
end

function assign_flowrefs!(c::MeasCollection)
    for analyzer in c.MoleAnalyzer
        assign_flowref!(analyzer, c.MassFlowMeas, c.VolumeFlowMeas, c.VolumeDensMeas)
    end
    for analyzer in c.MassAnalyzer
        assign_flowref!(analyzer, c.MassFlowMeas, c.VolumeFlowMeas, c.VolumeDensMeas)
    end
end

function calcvalues!(c::MeasCollection, thermo::ThermoModel)
    map(anylzr -> calcvalues!(anylzr, c, thermo), c.MoleAnalyzer)
    map(anylzr -> calcvalues!(anylzr, c, thermo), c.MassAnalyzer)
    return c 
end

function calcvalues!(analyzer::AbstractAnalyzer, c::MeasCollection, thermo::ThermoModel)
    flowref = analyzer.flowref

    if flowref.type <: VolumeFlowMeas
        calcvalues!(analyzer, c.VolumeFlowMeas[flowref], thermo)

    elseif flowref.type <: MassFlowMeas 
        calcvalues!(analyzer, c.MassFlowMeas[flowref], thermo)

    elseif flowref.type <: VolumeDensMeas 
        calcvalues!(analyzer, c.VolumeDensMeas[flowref], thermo)

    else 
        error("No mole flow calculation case for measurement type: $(flowref.type)")
    end
    return analyzer 
end



function loglik(x::AbstractVector, c::MeasCollection)
    return sum(fn-> loglik(x, c[fn]), fieldnames(MeasCollection))
end


function gettags(x::AbstractMeasCollection) 
    rawtags = TagInfo[]
    addtags!(rawtags, x)
    return validate_tags(rawtags)
end

function validate_tags(rawtags)
    sort!(rawtags, by=getname)
    validated_tags = [rawtags[begin]]

    for tag in rawtags[(begin+1):end]
        #Ensure the units are the same if there is a name collision
        if getname(tag) == getname(validated_tags[end])
            if getunits(tag) != getunits(validated_tags[end])
                error("Unit mismatch between $(tag) and $(validated_tags[end])")
            end
        else
            push!(validated_tags, tag)
        end
    end

    return validated_tags
end

function addtags!(tags::AbstractVector{TagInfo}, c::M) where M <: AbstractMeasCollection
    for fn in fieldnames(M)
        addtags!(tags, c[fn])
    end
    return tags
end

function addtags!(tags::AbstractVector{TagInfo}, vm::AbstractVector{<:AbstractMeas})
    for m in vm
        addtags!(tags, m)
    end
    return tags
end

addtags!(tags::AbstractVector{TagInfo}, m::AbstractMeas) = addtags!(tags, m.tag)
addtags!(tags::AbstractVector{TagInfo}, m::MoleBalance) = tags

function addtags!(tags::AbstractVector{TagInfo}, m::VolumeFlowMeas)
    addtags!(tags, m.statetag)
    addtags!(tags, m.tag)
    return tags 
end

function addtags!(tags::AbstractVector{TagInfo}, m::VolumeDensMeas)
    addtags!(tags, m.denstag)
    addtags!(tags, m.tag)
    return tags 
end

function addtags!(tags::AbstractVector{TagInfo}, newtags::AbstractVector{TagInfo})
    for newtag in newtags
        addtags!(tags, newtag)
    end
    return tags
end

function addtags!(tags::AbstractVector{TagInfo}, newtag::TagInfo)
    if hasname(newtag)
        push!(tags, newtag)
    end
    return tags
end
