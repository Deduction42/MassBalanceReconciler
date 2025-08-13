include("_AbstractMeas.jl")

abstract type AbstractFlowMeas{S, T} <: AbstractSingleMeas{S, T} end

#==========================================================================================================
Mass flow rates
==========================================================================================================#
@kwdef mutable struct MassFlowMeas{S, T, N} <: AbstractFlowMeas{S, T}
    id          :: Symbol
    tag         :: TagInfo
    value       :: T
    molarmass   :: Species{S, Float64, N}
    stream      :: StreamRef{S, N}
    stdev       :: Float64
end
MassFlowMeas{S, T}(x...) where {S,T}  = MassFlowMeas{S, T, length(S)}(x...)
MassFlowMeas{S, T}(;kw...) where {S,T}  = MassFlowMeas{S, T, length(S)}(;kw...)

function MassFlowMeas(measinfo::MeasInfo, streams::Dict{Symbol, <:StreamRef{S}}, thermo::ThermoModel) where S
    if length(measinfo.tags) != 1
        error("Measurement Type: MassFlowMeas only supports 1 tag, measurement id '$(measinfo.id)' contains $(length(measinfo.tags))")
    end

    tag = first(values(measinfo.tags))
    dimension(u"kg/s") == dimension(tag) || error("Tag $(tag) has units that are incompatible with mass flow rates")
    
    return MassFlowMeas{S, Float64}(
        id        = measinfo.id,
        tag       = tag,
        value     = NaN,
        stdev     = ustrip(dimension(u"kg/s"), measinfo.stdev),
        stream    = streams[measinfo.stream],
        molarmass = molar_weights(thermo)
    )
end

function prediction(x::AbstractVector, m::MassFlowMeas)
    stream = speciesvec(x, m.stream)
    return dot(m.molarmass[:], stream[:])
end


#==========================================================================================================
Mass flow rates
==========================================================================================================#
@kwdef mutable struct VolumeDensMeas{S, T, N} <: AbstractFlowMeas{S, T}
    id          :: Symbol
    tag         :: TagInfo
    value       :: T
    denstag     :: TagInfo
    densval     :: T
    molarmass   :: Species{S, Float64, N}
    stream      :: StreamRef{S, N}
    stdev       :: Float64
end
VolumeDensMeas{S, T}(x...) where {S,T}   = VolumeDensMeas{S, T, length(S)}(x...)
VolumeDensMeas{S, T}(;kw...) where {S,T} = VolumeDensMeas{S, T, length(S)}(;kw...)

function VolumeDensMeas(measinfo::MeasInfo, streams::Dict{Symbol, <:StreamRef{S}}, thermo::ThermoModel) where S
    if length(measinfo.tags) != 2
        error("Measurement Type: VolDensityMeas only supports 2 tags, measurement id '$(measinfo.id)' contains $(length(measinfo.tags))")
    end

    tags = measinfo.tags
    dimension(u"m^3/s")  == dimension(tags[:volflow]) || error("Tag $(tags[:volflow]) has units that are incompatible with volumetric flow rates")
    dimension(u"kg/m^3") == dimension(tags[:density]) || error("Tag $(tags[:density]) has units that are incompatible with density")
    
    return VolumeDensMeas{S, Float64}(
        id        = measinfo.id,
        tag       = tags[:volflow],
        value     = NaN,
        denstag   = tags[:density],
        densval   = 1.25, #Density of air at SATP as an initial value
        molarmass = molar_weights(thermo),
        stream    = streams[measinfo.stream],
        stdev     = ustrip(dimension(u"m^3/s"), measinfo.stdev)
    )
end

function readvalues!(m::VolumeDensMeas{S,T}, d::Dict, ind::Integer...) where {S,T}
    m.value = getvalue(d, m.tag, ind...)
    m.densval = getvalue(d, m.denstag, ind...)
    return m
end

function setvalues(m::M, v::T, d::T) where {T, L, M<:VolumeDensMeas{L}}
    function getfn(fn::Symbol) 
        if (fn == :value)
            return v
        elseif (fn == :densval)
            return d
        else 
            return getproperty(m, fn)
        end
    end 
    return basetype(M){L,T}(map(getfn, fieldnames(M))...)
end

function prediction(x::AbstractVector, m::VolumeDensMeas)
    stream = speciesvec(x, m.stream)
    return dot(m.molarmass[:], stream[:])/m.densval
end


#==========================================================================================================
Volumetric flow rates
==========================================================================================================#
@kwdef struct TempPress{E} <: FieldVector{2,E}
    T :: E 
    P :: E 
end

@kwdef mutable struct VolumeFlowMeas{S, T, N} <: AbstractFlowMeas{S, T}
    id       :: Symbol
    tag      :: TagInfo
    value    :: T
    statetag :: TempPress{TagInfo}
    stateval :: TempPress{T}
    molarvol :: Species{S, Float64, N}
    stream   :: StreamRef{S, N}
    stdev    :: Float64
end
VolumeFlowMeas{S, T}(x...) where {S,T} = VolumeFlowMeas{S, T, length(S)}(x...)
VolumeFlowMeas{S, T}(;kw...) where {S,T} = VolumeFlowMeas{S, T, length(S)}(;kw...)

function VolumeFlowMeas(measinfo::MeasInfo, streams::Dict{Symbol, <:StreamRef{S}}, thermo::ThermoModel) where S
    if length(measinfo.tags) != 3
        error("Measurement Type: VolumeFlowMeas only supports 3 tags, measurement id '$(measinfo.id)' contains $(length(measinfo.tags))")
    end
    N = length(S)
    stream = streams[measinfo.stream]
    tags = measinfo.tags
    Ttag = tags[:T]
    Ptag = tags[:P]

    #Assign default (SI) values to T,P if they are quantities, if they are tags, they will be overwritten
    Tval = hasname(Ttag) ? 298.15  : ustrip(u"K", getvalue(Ttag)*getunits(Ttag))
    Pval = hasname(Ptag) ? 101.3e3 : ustrip(u"Pa", getvalue(Ptag)*getunits(Ptag))

    #Check the dimension of the tag info 
    dimension(u"m^3/s") == dimension(tags[:V]) || error("Tag $(tags[:V]) has units that are incompatible with volumetric flow rates")

    thermostate = ThermoState{S, Float64}(
        model=thermo, 
        T=Tval, 
        P=Pval, 
        n=Species{S}(ones(N)./N),
        phase=stream.phase
    )

    return VolumeFlowMeas{S, Float64}(
        id       = measinfo.id,
        tag      = tags[:V], 
        value    = NaN,
        statetag = TempPress{TagInfo}(T=Ttag, P=Ptag),
        stateval = TempPress{Float64}(T=Tval, P=Pval),
        stdev    = ustrip(dimension(u"m^3/s"), measinfo.stdev),
        stream   = stream,
        molarvol = molar_volumes(thermostate)
    )
end

function readvalues!(m::VolumeFlowMeas{S,T}, d::Dict, ind::Integer...) where {S,T}
    m.stateval = TempPress(
        P = getvalue(d, m.statetag.P, ind...),
        T = getvalue(d, m.statetag.T, ind...)
    )
    m.value = getvalue(d, m.tag, ind...)
    return m
end

function setvalues(m::M, v::T, state::TempPress{T}) where {T, L, M<:VolumeFlowMeas{L}}
    function getfn(fn::Symbol) 
        if (fn == :value) 
            return v 
        elseif (fn == :stateval)
            return state
        else 
            return getproperty(m, fn)
        end
    end 
    return basetype(M){L,T}(map(getfn, fieldnames(M))...)
end

function prediction(x::AbstractVector, m::VolumeFlowMeas)
    stream = speciesvec(x, m.stream)
    return dot(m.molarvol[:], stream)
end

function updatethermo!(m::VolumeFlowMeas{S}, statevec::AbstractVector{<:Real}, thermo::ThermoModel) where S
    tp_state = m.stateval
    thermostate = ThermoState(
        model = thermo, 
        T = tp_state.T, 
        P = tp_state.P, 
        n = Species{S}(statevec[m.stream] .+ MOL_ε),
        phase = m.stream.phase
    )
    m.molarvol = molar_volumes(thermostate)
    return m
end


#==========================================================================================================
Mole/Mass flow rate calculations
==========================================================================================================#
function calculate_mole_flow(meas::VolumeFlowMeas, moles::Species, thermo::ThermoModel)
    molefracs = fractions(moles)
    thermostate = ThermoState(
        model = thermo, 
        T = meas.stateval[:T], 
        P = meas.stateval[:P], 
        n = molefracs,
        phase = meas.stream.phase
    )
    Mv = dot(molar_volumes(thermostate), molefracs)
    return meas.value/Mv
end

function calculate_mole_flow(meas::MassFlowMeas, molefracs::Species, thermo::ThermoModel)
    Mw = dot(meas.molarmass, fractions(molefracs))
    return meas.value/Mw
end

function calculate_mole_flow(meas::VolumeDensMeas, molefracs::Species, thermo::ThermoModel)
    Mw = dot(meas.molarmass, fractions(molefracs))
    return (meas.value*meas.densval)/Mw
end

function calculate_mass_flow(meas::VolumeFlowMeas, molefracs::Species, thermo::ThermoModel)
    Mw = dot(molar_weights(thermo), fractions(molefracs))
    return Mw*calculate_mole_flow(meas, molefracs, thermo)
end

function calculate_mass_flow(meas::MassFlowMeas, molefracs::Species, thermo::ThermoModel)
    return meas.value
end

function calculate_mass_flow(meas::VolumeDensMeas, molefracs::Species, thermo::ThermoModel)
    return meas.value*meas.densval
end