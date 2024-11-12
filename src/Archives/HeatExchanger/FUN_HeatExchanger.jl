
function predict_state!(hx1::AbstractHeatExchanger, hx0::AbstractHeatExchanger, ::Type{HX_Type}) where HX_Type <: HeatExchanger
    specs = hx0.specs

    function predict_temperature(stream, Q)
        if Q > 0 
            return stream.temperature + Q/(get_Cp(stream)*sum(stream.mass_flows))
        else
            return stream.temperature
        end
    end

    function predict_pressure(stream, A_Cf)
        return stream.pressure - pressure_drop(stream, A_Cf)
    end

    #Assign output identical to input
    hx1.heating_outlet.mass_flows = hx0.heating_inlet.mass_flows
    hx1.cooling_outlet.mass_flows = hx0.cooling_inlet.mass_flows
    hx1.R_Fouling = hx0.R_Fouling

    #Predict pressure
    hx1.heating_outlet.pressure = predict_pressure(hx0.heating_inlet, specs.hot_flow_A_Cf)
    hx1.cooling_outlet.pressure = predict_pressure(hx0.cooling_inlet, specs.cold_flow_A_Cf)

    #Temperature prediction (including zero flow cases)
    isRunning = (sum(hx0.heating_inlet.mass_flows)>0) && (sum(hx0.cooling_inlet.mass_flows)>0)
    α = exp( -MOD_ENV.ΔT*specs.UA/specs.Cp )

    if isRunning
        #Predict duty using current LMTD
        UA = 1/( 1/specs.UA + hx0.R_fouling )
        Q  = UA * hx0.LMTD

        #Predict temperature
        hx1.heating_outlet.temperature = predict_temperature(hx0.heating_inlet, -Q)
        hx1.cooling_outlet.temperature = predict_temperature(hx0.cooling_inlet, Q)

        #Heat exchanger Cp determines the dynamics of the equivalent LMTD
        hx1.LMTD = α*hx0.LMTD + (1-α)*log_mean_temperature_difference(hx0)        

    else #LMTD decays to nothing if not running
        hx1.LMTD = α*hx0.LMTD
    end

end

function log_mean_temperature_difference(hx::AbstractHeatExchanger)
    if hx.specs.is_counter_flow
        ΔT1  = hx.heating_inlet.temperature  - hx.cooling_outlet.temperature
        ΔT2  = hx.heating_outlet.temperature - hx.cooling_inlet.temperature
    else
        ΔT1  = hx.heating_inlet.temperature  - hx.cooling_inlet.temperature
        ΔT2  = hx.heating_outlet.temperature - hx.cooling_outlet.temperature        
    end

    return log_mean_temperature_difference(ΔT1, ΔT2)
end

function log_mean_temperature_difference(ΔT1::Real, ΔT2::Real)
    if sign(ΔT1) != sign(ΔT2)
        return 0 #LMTD requires the same sign, return mean
    end

    LMTD = (ΔT1-ΔT2) / log(ΔT1/ΔT2)
    if isnan(LMTD)
        return (ΔT1 + ΔT2)/2 #Just return mean if it doesn't work
    else
        return LMTD
    end
end
