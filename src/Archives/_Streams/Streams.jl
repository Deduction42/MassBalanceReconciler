
# ==============================================================================
# Process stream labels
# ==============================================================================
# We add a "label" to streams to prevent the proliferation of types, but it still allows for multiple dispatch
# Multiple dispatch allows for differnet methods like "get density" on different stream types
abstract type Substance end

#CustomSpec means you define all basic thermodynamic properties yourself
abstract type CustomSpec <: Substance end 

abstract type AbstractLiquid <: Substance end
abstract type AbstractGas <: Substance end
abstract type AbstractVLE <: Substance end 
abstract type AbstractSlurry <: Substance end
abstract type AbstractSolid  <: Substance end

abstract type IdealVLE    <: AbstractVLE end
abstract type IdealLiquid <: AbstractLiquid end
abstract type IdealGas    <: AbstractGas end
abstract type IdealSlurry <: AbstractSlurry end



# ==============================================================================
# Process stream object
# ==============================================================================
#Custom properties that can be specified by shortcut
Base.@kwdef struct StreamSpecs <: AbstractSpecs
    resistance :: Float32 = 0.0
    components :: Vector{Symbol} = CHEM.Names
    density :: Float32 = NaN
    Cp :: Float32 = NaN
    viscosity :: Float32 = NaN
end

mutable struct ProcessStream{T} <: AbstractConnector{T}
    @inherit_asset_fields(GenericAsset{T})
    subclass :: DataType
    temperature :: T
    pressure  :: T
    components :: Vector{T}
    specs :: StreamSpecs
end

#Default values of streams, boundaries and specs ------------------------------------------------------
function default_value(::Type{asset}, fn::Symbol) where asset <: ProcessStream{<:Number}
    ϵ = 1e-32; #Very small number as many things break down at 0
    T = eltype(asset)

    if fn == :temperature
        return T(298.15)
    elseif fn == :pressure 
        return T(101.0)
    elseif fn == :components
        return fill(T(ϵ), length(CHEM.Names))
    else
        error("Asset type \"$(asset)\" has no default value for field \"$(fn)\"")
    end
end

function default_value(::Type{asset}, fn::Symbol) where asset <: ProcessStream{Bounds}
    ϵ = 1e-32; #Very small number as many things break down at 0

    if fn == :temperature #Kelvin
        return Bounds(ϵ,1e7)
    elseif fn == :pressure #kPa
        return Bounds(ϵ,1e7) 
    elseif fn == :components #kg or kg/s
        return fill(Bounds(ϵ,1e9), length(CHEM.Names))
    else
        error("Asset type \"$(asset)\" has no default value for field \"$(fn)\"")
    end
end

function build_specs(::Type{specType}, D::Nothing) where specType <: Union{ProcessStream, StreamSpecs}
    return StreamSpecs()
end

function build_specs(::Type{specType}, D::Dict{String}) where specType <: Union{ProcessStream, StreamSpecs}
    defaultSpec = StreamSpecs()
    
    function fallback_get(fn::Symbol) 
        return if fn == :components
            Symbol.(get(D, string(fn), getproperty(defaultSpec, fn)))
        else
            get(D, string(fn), getproperty(defaultSpec, fn))
        end
    end

    return StreamSpecs( (fallback_get(fn) for fn in fieldnames(StreamSpecs))... )
end

#Easy construction ------------------------------------------------------
function ProcessStream(x; substance=Substance)
    n = length(CHEM.Names)
    T = typeof(x)
    return ProcessStream{T}(:NA, substance, x, x, fill(x,n), StreamSpecs())
end

function fill_asset(x, oldStream::ProcessStream)
    newStream = ProcessStream(x)
    for fn in passthrough_fields(typeof(oldStream))
        newStream[fn] = oldStream[fn]
    end
    return newStream
end




#Lists the numeric indices of the available components in a stream 
#Some streams don't have certain components (so this forces zero flows for some components)
#Fallback method which returns all indices of the components
function stream_components(s::ProcessStream)
    return stream_components(s.subclass)
end

function stream_components(::Type{T}) where T <: Substance
    return collect( eachindex(CHEM.Names) )
end

function mix(s1::ProcessStream, s2::ProcessStream; CombinedSubtance=nothing)
    s3 = deepcopy(s1)
    s3.subclass = something(CombinedSubstance, s1.subclass)

    return mix!(s3, s1, s2)
end

function mix!(s3::ProcessStream, s1::ProcessStream, s2::ProcessStream)
    s3.pressure = min(s1.pressure, s2.pressure)
    s3.components .= s1.components .+ s2.components
    s3.temperature = mix_temperature(s1, s2)

    return s3
end


function mix_temperature(s1::ProcessStream, s2::ProcessStream)
    mCp(x) = get_Cp(x)*total_mass(x)
    (mCp1, mCp2) = (mCp(s1), mCp(s2))
    return (mCp1*s1.temperature + mCp2*s2.temperature)/(mCp1 + mCp2)
end

# ==============================================================================
# Special case state handling for streams
# ==============================================================================
function key_states(::Type{T}) where T <: ProcessStream
    return (:temperature, :pressure, :components)
end


function set_state!(state::ProcessStream, stateVec::AbstractVector, ind::ProcessStream)
    state.temperature = stateVec[ind.temperature]
    state.pressure    = stateVec[ind.pressure]

    for (ii,v) in enumerate(ind.components)
        if v != 0
            state.components[ii] = stateVec[v]
        else
            state.components[ii] = 0
        end
    end

    return state
end

function set_state_deep!(state::ProcessStream, stateVec::AbstractVector, ind::ProcessStream)
    return set_state!(state, stateVec, ind)
end

function get_state!(stateVec::AbstractVector, state::ProcessStream, ind::ProcessStream)
    stateVec[ind.temperature] = state.temperature
    stateVec[ind.pressure] = state.pressure

    for (ii,v) in enumerate(ind.components)
        if v != 0
            stateVec[v] = state.components[ii]
        end
    end

    return stateVec
end


function fill_indexer!(lastInd::Index_Type, indexer::ProcessStream; CompIndexer=COMPONENT_INDEXER)
    indexer.temperature = iterate_indexer!(lastInd)
    indexer.pressure    = iterate_indexer!(lastInd)

    #Fill in the active components (indexer for inactive is 0 by default)
    indexer.components .= 0
    for ci in (CompIndexer[c] for c in indexer.specs.components)
        indexer.components[ci] = iterate_indexer!(lastInd)
    end

    return indexer
end


# ==============================================================================
# Key functions for calculated states for streams
# ==============================================================================
function total_mass(X::ProcessStream)
    return sum(X.components)
end

fracvec(x::Vector{<:Number}) = x ./ nonzero(sum(x))

mass_fractions(X::ProcessStream) = fracvec(X.components)
mole_fractions(X::ProcessStream) = fracvec(X.components ./ CHEM.MolWeight)

mass2molefrac(x::Vector{<:Number}) = fracvec(x ./ CHEM.MolWeight)
mole2massfrac(x::Vector{<:Number}) = fracvec(x .* CHEM.MolWeight)

get_density(X::ProcessStream) = 1 / nonzero( get_specific_volume(X) )

get_specific_volume(X::ProcessStream) = stream_dispatch(get_specific_volume, X)
get_specific_volume(::Type{T}, X::ProcessStream) where T<:Substance = dot(mass_fractions(X), component_specific_volume(T,X))

component_specific_volume(X::ProcessStream) = stream_dispatch(component_specific_volume, X)



function pressure_drop(X::ProcessStream, R)
    ρ = get_density(X)
    m = sum(X.components)

    return (R*m)^2 / (2*ρ)
end

function component_specific_volume(::Type{T}, X::ProcessStream) where {T <: IdealGas}
    R = CHEM.R
    n = 1 ./ CHEM.MolWeight

    return n.*R.*(X.temperature) ./ X.pressure
end

function component_specific_volume(::Type{T}, X::ProcessStream) where {T <: IdealLiquid}
    return copy(CHEM.LiqVol)
end

function component_specific_volume(::Type{T}, X::ProcessStream) where {T <: IdealSlurry}
    specVol = copy(CHEM.LiqVol)

    for ii in 1:length(specVol)
        if (MOD_ENV.ComponentMeltPoint[ii]) > X.temperature
            specVol[ii] = MOD_ENV.ComponentSolVol[ii]
        end
    end

    return specVol
end


#=
# ==============================================================================
# Special case, water stream
# ====================================c==========================================
abstract type ProcessStream_Water <: AbstractStreamLabel end

function stream_components(::Type{T}; ci=COMPONENT_INDEXER) where T <: ProcessStream_Water
    return [ci.H2O]
end

function get_density(X::ProcessStream, ::Type{T}) where T <: ProcessStream_Water
    return 1000 #kg/m3
end

function get_Cp(X::ProcessStream, ::Type{T}) where T <: ProcessStream_Water
    return 4.2  #kJ/kg
end

function get_specific_volumes(X::ProcessStream, ::Type{T}) where T <: ProcessStream_Water
    volumes = deepcopy(X.components)

    for fn in propertynames(Volumes)
        if fn == :H2O
            volumes[fn] = 0.001
        else
            volumes[fn] = 0
        end
    end

    return volumes
end

# ==============================================================================
# Special case, ideal gas
# ==============================================================================
abstract type ProcessStream_IdealGas <: AbstractStreamLabel end


function get_specific_volumes(X::ProcessStream, ::Type{T}) where T <: ProcessStream_IdealGas
    volumes  = (8.314*(X.temperature+273)/X.pressure) ./ StreamComponentRef.molecular_weights
    replace!(volumes, NaN=>0)

    return volumes
end

function get_density(X::ProcessStream, ::Type{T}) where T <: ProcessStream_IdealGas
    specVol = get_specific_volumes(X)
    ρ = sum(X.components) / sum( specVol .* X.components )
    return ρ
end
=#

#= Testing
Source = ProcessStream{Float64}(Default=0, StreamType=ProcessStream_IdealGas)
Source.id = :S1

Target = ProcessStream{Float64}(StreamType=ProcessStream_IdealGas)
Target.id = :S2

streamInd = fill_indexer_deep!(Index_Type(0), Source)
update_states_deep!(Source, [0.5, 1.0, 1.5, 205, 300], streamInd )
Target = copy_into_deep!(Target, Source)

=#
