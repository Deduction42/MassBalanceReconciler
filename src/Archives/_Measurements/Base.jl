abstract type AbstractMeasurement{T} <: AbstractAsset{T} end

struct ConversionToSI
    mult :: Float64
    add  :: Float64
end
ConversionToSI(x::AbstractVector) = ConversionToSI(x[1], x[2])
ConversionToSI(x::Number) = ConversionToSI(x, 0)

mutable struct BaseMeasurement{T} <: AbstractMeasurement{T}
    id   :: Symbol
    asset_id :: Symbol
    value :: T
    stdev :: Float64
    to_SI :: ConversionToSI
end

function default_value(::Type{asset}, fn::Symbol) where asset <: AbstractMeasurement{<:Number}
    T = eltype(asset)
    if fn == :value
        return T(0)
    elseif fn == :to_SI 
        return ConversionToSI(1,0)
    else
        error("Asset type \"$(asset)\" has no default value for field \"$(fn)\"")
    end
end

get_si_measurement(m::AbstractMeasurement)    = convert_to_si(m.value, m.to_SI)
convert_to_si(x::Number, c::ConversionToSI)   = c.mult*x + c.add
convert_from_si(x::Number, c::ConversionToSI) = (x - c.add)/c.mult


function predict_measurements!(res::AbstractVector, measGroups::NamedTuple, BOM::NamedTuple)
    f(meas::AbstractMeasurement) = predict_measurement(meas, BOM[meas.asset_id])
    return collect_measgroup_results!(f, res, measGroups)
end

function get_si_measurements!(res::AbstractVector, measGroups)
    return collect_measgroup_results!(get_si_measurement, res, measGroups)
end

function get_meas_variances!(res::AbstractVector, measGroups)
    f(meas::AbstractMeasurement) = meas.stdev^2
    return collect_measgroup_results!(f, res, measGroups)
end

#Collect results from measurement groups without allocation, the result vector is the same as the sum of the measgroup lengths
function collect_measgroup_results!(f::Function, res::AbstractVector, measGroups::NamedTuple)
    ii = 1
    for measVec in values(measGroups)
        for meas in measVec
            res[ii] = f(meas)
            ii += 1
        end
    end
    return res
end





#=
function predict_measurements(measGroups::NamedTuple, BOM::NamedTuple)
    #Generate a prediction vector of all measurements (grouped by measrement type)
    vecGen = ( predict_measurements(measVec, BOM) for measVec in values(measGroups) )
    return vcat( vecGen... )
end    

function predict_measurements!(measList::Vector{<:AbstractMeasurement}, BOM)
    #Generate a prediction vector of each measurement
    return [ predict_measurement(meas, BOM[meas.asset_id]) for meas in measList ]
end
=#





#= STRATEGY =================================================================================
to_SI applies to measurements coming IN (state = meas.value*meas.scale_factor)
An object "measurementBOM" will be a NamedTuple of measruement objects
Each entry of the measruementBOM will contain a vector of a measurement type
This allows for a minimal type signature and an organization of the measurement vectors
=#

function collect_by_type(measList)
    measTypes  = unique(typeof.(measList))
    measGroups = (; (Symbol(T)=>T[] for T in measTypes)... )

    for meas in values(measList)
        push!( measGroups[Symbol(typeof(meas))], meas)
    end

    return measGroups
end



