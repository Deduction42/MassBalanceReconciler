function predict_state!(sep1::ObjType, sep0::ObjType, ::Type{AsType}) where {ObjType<:AbstractSeparatorVLE, AsType<:SeparatorVLE}
    sep1.inlet = deepcopy(sep0.inlet)

    sep1.vapor_frac = copy(sep0.vapor_frac) #Random walk vapor prediction
    mass_balance_prediction!(sep1, sep0, SeparatorVLE)
    energy_balance_prediction!(sep1, sep0, SeparatorVLE)

    return nothing
end


#Predicts output flows based on inputs
function mass_balance_prediction!(sep1::ObjType, sep0::ObjType, ::Type{AsType}) where {ObjType<:AbstractSeparatorVLE, AsType<:SeparatorVLE}
    Δt = MOD_ENV.ΔT

    #Fraction of output leaving in gas form (use sep1 which may contain predictions for this)
    γ = sep1.vapor_frac

    #Average mass over time period
    Σm  = sum(sep0.contents.components) + 0.5*sep0.accumulation*Δt
    totalOutFlow = sum(sep0.inlet.contents) - sep0.accumulation*Δt

    #Integration of contents to find step change (dynamic simulation due to concentrations)
    τ = Σm ./ totalOutFlow
    α = exp.(-Δt./τ)
    Δcontents = (α.-1).*sep0.contents.components + τ.*(α.-1).*sep0.inlet.components

    #Update contents
    sep1.contents.components = sep0.contents.components + Δcontents

    #Predict output flow rates
    outFlows = sep0.inlet.components .- (Δcontents./Δt)
    sep1.gas_outlet.contents .= outFlows .* γ
    sep1.liquid_outlet.contents .= outFlows .* (1 .- γ) 

    return nothing
end

#Predicts duty (incoming energy or outgoing energy) after mass balance has been completed
function energy_balance_prediction!(sep1::ObjType, sep0::ObjType, ::Type{AsType}) where {ObjType<:AbstractSeparatorVLE, AsType<:SeparatorVLE}
    Δt  = MOD_ENV.ΔT

    sep1.gas_outlet.temperature    = sep1.contents.temperature
    sep1.gas_outlet.pressure       = sep1.conents.pressure
    sep1.liquid_outlet.temperature = sep1.contents.temperature
    sep1.liquid_outlet.pressure    = sep1.conents.pressure

    Qin  = get_mCp(sep0.inlet).*(sep1.contents.temperature - sep0.inlet.temperature)
    Qin += dot(sep1.gas_outlet.components, component_evap_heat(sep1.gas_outlet))

    sep1.heat_input = Qin

    return nothing
end
