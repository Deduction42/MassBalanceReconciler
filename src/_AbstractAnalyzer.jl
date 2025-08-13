include("_AbstractFlowMeas.jl")

abstract type AbstractAnalyzer{S, T} <: AbstractMultiMeas{S, T} end

prediction(x::AbstractVector, m::AbstractAnalyzer) = speciesvec(x, m.stream)
getvalue(m::AbstractAnalyzer) = speciesvec(m.moleflow)


function assign_flowref!(m::AbstractAnalyzer, measvecs::AbstractVector{<:AbstractFlowMeas}...)
    m.flowref = MeasReference(m.flowref.tag, measvecs...)
    return m 
end


#==========================================================================================================
Molar Analysis
==========================================================================================================#
@kwdef mutable struct MoleAnalyzer{S, T, N} <: AbstractAnalyzer{S, T}
    id        :: Symbol
    tag       :: Species{S, TagInfo, N}
    molefrac  :: Species{S, T, N}
    moleflow  :: Species{S, T, N}
    stream    :: StreamRef{S, N}
    stdev     :: Species{S, Float64, N}
    flowref   :: MeasReference
end
MoleAnalyzer{S, T}(x...) where {S,T} = MoleAnalyzer{S, T, length(S)}(x...)
MoleAnalyzer{S, T}(;kw...) where {S,T} = MoleAnalyzer{S, T, length(S)}(;kw...)

function MoleAnalyzer(measinfo::MeasInfo, streams::Dict{Symbol, <:StreamRef{S}}, thermo::ThermoModel) where S
    measid = measinfo.id
    N = length(S)

    if length(measinfo.tags) != (N+1)
        error("Measurement Type: MassFlowMeas only supports $(N+1) tags, measurement id '$(measid)' contains $(length(measinfo.tags))")
    end

    #Initialize a junk value for now, will be replaced later
    flowref = MeasReference(
        tag  = measinfo.tags[:TOTAL_FLOW],
        type = AbstractFlowMeas,
        ind  = 0
    )

    moleflow_σ = zero(SVector{N,Float64}) .+ ustrip(dimension(u"mol/s"), measinfo.stdev)

    return MoleAnalyzer{S, Float64}(
        id       = measid,
        tag      = Species{S, TagInfo, N}(measinfo.tags),
        molefrac = zero(Species{S, Float64, N}),
        moleflow = zero(Species{S, Float64, N}),
        stream   = streams[measinfo.stream],
        stdev    = Species{S, Float64, N}(moleflow_σ),
        flowref  = flowref
    )
end

function readvalues!(m::MoleAnalyzer, d::Dict, ind::Integer...)
    m.molefrac = fractions(getvalue(d, m.tag, ind...))
    return m 
end

function setvalues(m::M, mfrac::Species{L,T}, mflow::Species{L,T}) where {T, L, M<:MoleAnalyzer{L}}
    function getfn(fn::Symbol) 
        if (fn == :molefrac)
            return mfrac
        elseif (fn == :moleflow)
            return mflow
        else 
            return getproperty(m, fn)
        end
    end
    return basetype(M){L,T}(map(getfn, fieldnames(M))...)
end

function calcvalues!(analyzer::MoleAnalyzer, flow::AbstractFlowMeas, thermo::ThermoModel)
    totalflow = calculate_mole_flow(flow, analyzer.molefrac, thermo)
    analyzer.moleflow = totalflow.*(analyzer.molefrac)
    return analyzer 
end


#==========================================================================================================
Mass Analysis
==========================================================================================================#
@kwdef mutable struct MassAnalyzer{S, T, N} <: AbstractAnalyzer{S, T}
    id        :: Symbol
    tag       :: Species{S, TagInfo, N}
    massfrac  :: Species{S, T, N}
    moleflow  :: Species{S, T, N}
    stream    :: StreamRef{S, N}
    stdev     :: Species{S, Float64, N}
    flowref   :: MeasReference
end
MassAnalyzer{S, T}(x...) where {S,T} = MassAnalyzer{S, T, length(S)}(x...)
MassAnalyzer{S, T}(;kw...) where {S,T} = MassAnalyzer{S, T, length(S)}(;kw...)

function MassAnalyzer(measinfo::MeasInfo, streams::Dict{Symbol, <:StreamRef{S}}, thermo::ThermoModel) where S
    measid = measinfo.id
    N = length(S)

    if length(measinfo.tags) != (N+1)
        error("Measurement Type: MassFlowMeas only supports $(N+1) tags, measurement id '$(measid)' contains $(length(measinfo.tags))")
    end

    #Initialize a junk value for now, will be replaced later
    flowref = MeasReference(
        tag  = measinfo.tags[:TOTAL_FLOW],
        type = AbstractFlowMeas,
        ind  = 0
    )

    moleflow_σ = zero(SVector{N,Float64}) .+ ustrip(dimension(u"mol/s"), measinfo.stdev)

    return MassAnalyzer{S, Float64}(
        id       = measid,
        tag      = Species{S, TagInfo, N}(measinfo.tags),
        massfrac = zero(Species{S, Float64, N}),
        moleflow = zero(Species{S, Float64, N}),
        stream   = streams[measinfo.stream],
        stdev    = Species{S, Float64, N}(moleflow_σ),
        flowref  = flowref
    )
end


function readvalues!(m::MassAnalyzer, d::Dict, ind::Integer...)
    m.massfrac = fractions(getvalue(d, m.tag, ind...))
    return m 
end

function setvalues(m::M, mfrac::Species{L,T}, mflow::Species{L,T}) where {T, L, M<:MassAnalyzer{L}}
    function getfn(fn::Symbol) 
        if (fn == :massfrac)
            return mfrac
        elseif (fn == :moleflow)
            return mflow
        else 
            return getproperty(m, fn)
        end
    end
    return basetype(M){L,T}(map(getfn, fieldnames(M))...)
end

function calcvalues!(analyzer::MassAnalyzer{S}, flow::AbstractFlowMeas, thermo::ThermoModel) where S
    molefracs = Species{S}(fractions(analyzer.massfrac ./ molar_weights(thermo)))
    totalflow = calculate_mole_flow(flow, molefracs, thermo)
    analyzer.moleflow = totalflow.*(molefracs)
    return analyzer 
end




