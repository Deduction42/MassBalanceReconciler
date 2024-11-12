using QuadGK

# ====================================================================================
# Enthalpy of evaporation
# ====================================================================================
function stream_ΔH_vap(st::ProcessStream)
    sum( st.components .* component_ΔH_vap(st) )
end

function component_ΔH_vap(st::ProcessStream)
    n = length(st.components)
    ΔHv = zeros(eltype(st), n)

    for (ii,mc) in enumerate(st.components)
        if mc > 0
            #Reduced temperature correlation from Watson's equation
            Tr = min(1, st.temperature/CHEM.CriticalTemperature[ii])
            ΔHv[ii] = calc_ΔHv(Tr, CHEM.ΔH_vap[ii])
        end
    end

    return ΔHv
end

function calc_ΔH_vap(Tr::Number, params::Polynomial)
    C = params.Params
    return C[1]*(1-Tr)^(C[2] + C[3]*Tr + C[4]*Tr^2)
end

# ====================================================================================
# Enthalpy of fusion
# ====================================================================================
function stream_ΔH_fus(st::ProcessStream)
    sum( st.components .* component_ΔH_fus(st) )
end

function component_ΔH_fus(st::ProcessStream)
    n = length(st.components)
    ΔHf = zeros(eltype(st), n)

    for (ii,mc) in enumerate(st.components)
        if mc > 0
            #Reduced temperature correlation from Watson's equation
            Tr = min(1, st.temperature/CHEM.CriticalTemperature[ii])
            ΔHf[ii] = CHEM.ΔH_fus[ii]
        end
    end

    return ΔHf
end

# ====================================================================================
# Heat capacities
# ====================================================================================
function stream_mCp(X::ProcessStream)
    return sum( X.components .* component_mCp(X) )
end

function component_mCp(X::ProcessStream)
    return phase_dispatch(calc_mCp, X.subclass, X.components, X.temperature)
end

function component_∫mCp_dT(X::ProcessStream, ΔT::Pair{<:Number, <:Number})
    return phase_dispatch(calc_∫mCp_dT, X.subclass, X.components, Tuple(ΔT))
end

function phase_dispatch(callfunc::Function, phase::DataType, argList...)
    phaseList = (AbstractGas, AbstractLiquid, AbstractSolid, AbstractSlurry, AbstractVLE)

    if phase <: AbstractGas
        return callfunc(AbstractGas, argList...)
    elseif phase <: AbstractLiquid
        return callfunc(AbstractLiquid, argList...)
    elseif phase <: AbstractSolid
        return callfunc(AbstractSolid, argList...)
    elseif phase <: AbstractSlurry
        return callfunc(AbstractSlurry, argList...)
    elseif phase <: AbstractVLE
        return callfunc(AbstractVLE, argList...)   
    else
        println("WARNING: Stream type $(phase) doesn't have a phase classification, using generic fallback")
        return callfunc(Substance, argList...)     
    end
end


# ====================================================================================
# Calculate the heat capcacity (non-specific)
# ====================================================================================
SinglePhase = Union{AbstractGas, AbstractLiquid, AbstractSolid}

#Single phase --------------------------------------------------------------------------------------------------------
function calc_mCp(::Type{phase}, st::ProcessStream, temperature::Number) where phase <: SinglePhase
    return calc_mCp(phase, st.components, temperature)
end

function calc_mCp(::Type{phase}, massVec::Vector{<:Number}, temperature::Number) where phase <: SinglePhase
    mCp = zeros(eltype(massVec), length(massVec));
    for (ii, mc) in enumerate(massVec) 
        if mc > 0
            mCp[ii] = (mc/CHEM.MolWeight[ii]) * calc_Cp(phase, ii, temperature)
        end
    end
    return mCp #kJ/K
end

function calc_∫mCp_dT(::Type{phase}, st::ProcessStream, temperature::NTuple{2,<:Number}) where phase <: SinglePhase
    return calc_∫mCp_dT(phase, st.components, temperature)
end

function calc_∫mCp_dT(::Type{phase}, massVec::Vector{<:Number}, temperature::NTuple{2,<:Number}) where phase <: SinglePhase
    mCpΔT = zeros(eltype(massVec), length(massVec));
    for (ii, mc) in enumerate(massVec)
        if mc > 0
            mCpΔT[ii] = (mc/CHEM.MolWeight[ii]) * calc_∫Cp_dt(phase, ii, temperature)
        end
    end
    return mCpΔT #kJ
end


#Slurries ------------------------------------------------------------------------------------------------------------
function slurry_Cp(func_Cp, st::ProcessStream, temperature::Number)
    liqInd = (st.temperature .>= CHEM.ComponentMeltPoint)
    return +(
        func_Cp(AbstractLiquid, st.components.*liqInd, temperature),
        func_Cp(AbstractSolid,  st.components.*(.!liqInd), temperature)
    )
end

function calc_mCp(::Type{phase}, st::ProcessStream, temperature::Number) where phase <: AbstractSlurry
    return slurry_Cp(calc_mCp, st, temperature)
end

function calc_∫mCp_dT(::Type{phase}, st::ProcessStream, temperature::NTuple{2,<:Number}) where phase <: AbstractSlurry
    return slurry_Cp(calc_∫mCp_dT, st, temperature)
end

#Vapor Liquid Equilibrium -------------------------------------------------------------------------------------------
function vle_Cp(func_Cp, st::ProcessStream, temperature::Number)
    gasFrac = gas_fractions(st)
    return +(
        func_Cp(AbstractLiquid, st.components.*(1.0.-gasFrac), temperature),
        func_Cp(AbstractSolid,  st.components.*gasFrac, temperature)
    )
end

function calc_mCp(::Type{phase}, st::ProcessStream, temperature::Number) where phase <: AbstractVLE
    return vle_Cp(calc_mCp, st, temperature)
end

function calc_∫mCp_dT(::Type{phase}, st::ProcessStream, temperature::NTuple{2,<:Number}) where phase <: AbstractVLE
    return vle_Cp(calc_∫mCp_dT, st, temperature)
end

#Last-ditch fallback (phase agnostic) -------------------------------------------------------------------------------
function fallback_Cp(func_Cp, st::ProcessStream, temperature::Number)
    return +(
        func_Cp(AbstractGas, st.components./3, temperature),
        func_Cp(AbstractLiquid, st.components./3, temperature),
        func_Cp(AbstractSolid,  st.components./3, temperature)
    )
end

function calc_mCp(::Type{phase}, st::ProcessStream, temperature::Number) where phase <: Substance
    return fallback_Cp(calc_mCp, st, temperature)
end

function calc_∫mCp_dT(::Type{phase}, st::ProcessStream, temperature::NTuple{2,<:Number}) where phase <: Substance
    return fallback_Cp(calc_∫mCp_dT, st, temperature)
end


# ====================================================================================
# Calculate the specific heat capacity, dispatch on Cas/Liquid/Solid phases
# ====================================================================================
calc_Cp(::Type{phase}, ind::Integer, T::Number) where phase <: AbstractSolid  = poly_fx(T, CHEM.Cp_Sol[ind])
calc_Cp(::Type{phase}, ind::Integer, T::Number) where phase <: AbstractLiquid = poly_fx(T, CHEM.Cp_Liq[ind])
calc_Cp(::Type{phase}, ind::Integer, T::Number) where phase <: AbstractGas    = calc_CpGas(T, CHEM.Cp_Gas[ind])

calc_∫Cp_dt(::Type{phase}, ind::Integer, T::NTuple{2}) where phase <: AbstractSolid  = poly_∫fdx(T, CHEM.Cp_Sol[ind])
calc_∫Cp_dt(::Type{phase}, ind::Integer, T::NTuple{2}) where phase <: AbstractLiquid = poly_∫fdx(T, CHEM.Cp_Liq[ind])
calc_∫Cp_dt(::Type{phase}, ind::Integer, T::NTuple{2}) where phase <: AbstractGas    = calc_∫CpGas(T, CHEM.Cp_Gas[ind])


# Hyperbolic correlation for gasses --------------------------------------------------------------
function calc_CpGas(T::Number, params::Polynomial)
    C = params.Params
    C3Ti = C[3]/T
    C5Ti = C[5]/T
    return C[1] + C[2]*(C3Ti/sinh(C3Ti))^2 + C[4]*(C5Ti/cosh(C5Ti))^2
end

function calc_∫CpGas(T::NTuple{2,<:Number}, params::Polynomial)
    (vT, w) = gauss(6, T[1], T[2]) #Gauss quadrature points
    return sum( calc_CpGas.(vT, [params]) .* w )
end

# Solids have more or less a constant Cp unless it gets very cold 
# Cp = 3nR kJ/kmol-K where "n" is the number of atoms per mole
# http://vallance.chem.ox.ac.uk/pdfs/EinsteinDebye.pdf





