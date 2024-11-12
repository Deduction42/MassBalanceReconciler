abstract type WaterStream <: AbstractLiquid end

#Specific volume equation (it's essentially a pure liquid so it's ideal)
function get_specific_volume(::Type{T}, x::ProcessStream) where T <: WaterStream
    return 0.001
end