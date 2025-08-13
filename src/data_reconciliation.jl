include("_PlantState.jl")
#using Optim, LineSearches
using Optimization
import Zygote
import OSQP


#==========================================================================================================
Observer matrix construction
==========================================================================================================#
@kwdef struct ObserverView{V<:AbstractVector, M<:AbstractMatrix}
    stateind :: V 
    jacobian :: M 
end

function ObserverView(m::AbstractMeas, x::AbstractVector)
    as_vec(x::Number) = SVector(x)
    as_vec(x::AbstractVector) = x

    ind = state_indices(m)

    function abridged_observer(xi)
        spx = sparsevec(ind, xi)
        return as_vec(prediction(spx, m))
    end

    obsview = ObserverView(
        stateind = ind,
        jacobian = ForwardDiff.jacobian(abridged_observer, x[ind])
    )

    return obsview
end



function reconcile!(plant::PlantState, interval::TimeInterval, data::AbstractDict{K, TimeSeries{T}}) where {K<:AbstractString, T<:Real}
    get_statevec(plant::PlantState) = deepcopy(plant.statevec)
    get_statestd(plant::PlantState) = sqrt.(inv.(Vector(diag(plant.stateinv))))

    #Average the data over the intervals
    data_avg = si_time_averages(plant, interval, data)
    vt = data_avg[TIMESTAMP_KEY]

    #Initialize the samples and states
    samples = Dict{String,Float64}()
    t0 = interval[begin]
    states  = TimeSeries([TimeRecord(t0, get_statevec(plant))])
    stdevs  = TimeSeries([TimeRecord(t0, get_statestd(plant))])

    #Reconciliation for each timestamp
    for ii in eachindex(vt)
        for (k, v) in pairs(data_avg)
            samples[k] = v[ii]
        end
        t  = samples[TIMESTAMP_KEY]
        optimresults = reconcile!(plant, samples)

        @info "$(round(ii/length(vt)*100)) % complete"

        push!(states, TimeRecord(t, get_statevec(plant)))
        push!(stdevs, TimeRecord(t, get_statestd(plant)))
    end

    return PlantSeries(plant, states, stdevs)
end


"""
si_time_averages(plant::PlantState, interval::TimeInterval, data::AbstractDict{<:String, TimeSeries{T}}) where T

Use numerical integration to find averages on 'data', over the time interval 'interval' at the desired sampling rate 
defined in 'plant.clock.interval'. This produces a regularly sampled Dict{String, Vector{T}} where values are in SI units
"""
function si_time_averages(plant::PlantState, interval::TimeInterval, data::AbstractDict{<:String, TimeSeries{T0}}) where T0
    #Create a vector of sampled timestamps
    Δt = round(plant.clock.interval[])
    vt = interval[begin]:Δt:interval[end]

    #Obtain tags and initialize the averages
    T = promote_type(T0, Float64)
    data_avg = Dict{String, Vector{T}}()

    for taginfo in plant.tags #Calculate the averages and convert to SI units
        vec_avg  = values(average(data[taginfo.name], vt, order=0))
        vec_avg .= to_si_units.(vec_avg, taginfo.units)
        data_avg[taginfo.name] = vec_avg
    end

    #Timestamps representing the end of the averaging period
    data_avg[TIMESTAMP_KEY] = vt[(begin+1):end]

    #Zero measurement tag 
    data_avg[""] = fill(convert(T, 0.0), length(vt)-1)

    return data_avg
end

to_si_units(x, u::AbstractUnits)   = (u |> dimension(u))(x)
from_si_units(x, u::AbstractUnits) = (dimension(u) |> u)(x)


function reconcile!(plant::PlantState, data::AbstractDict{<:String, <:Real})
    readvalues!(plant, data)
    updatethermo!(plant)
    predict!(plant)

    optimresults = reconcile_plant!(plant)
    update_balance_errors!(plant)
    return optimresults
end

function reconcile_plant!(plant::PlantState; maxiter=3)
    measurements  = plant.measurements
    x   = copy(plant.statevec)
    y   = measurement_vector(measurements)
    N   = length(x)
    lb  = zeros(N)
    ub  = fill(1e9, N)

    #Replace NaNs with previous state estimate, dramatically increase varaince
    nan_ind = isnan.(y)
    σ² = variance_vector(measurements)
    σ²[nan_ind] .=  σ²[nan_ind] .* 100
    y[nan_ind]  .= prediction_vector(measurements, x)[nan_ind]

    #Closure to calculate the root-mean-square-error
    function calculate_rmse(x)
        updatethermo!(plant)
        ε² = abs2.(y .- prediction_vector(measurements, x))
        return sqrt(sum(ε²./σ²)/length(σ²))
    end
    
    rmse = Ref(calculate_rmse(x))

    for ii in 1:maxiter
        P⁻¹ = solver_setup!(plant)
        results = OSQP.solve!(plant.solver)
        new_rmse = calculate_rmse(results.x)

        if (new_rmse > rmse[])
            return results
        else
            x .= clamp.(results.x, lb, ub)
            plant.statevec .= x 
            
            if ii == maxiter
                plant.stateinv .= P⁻¹
                return results
            else 
                rmse[] = new_rmse
            end
        end
    end

    error("Uncaught branch in function")
end

function solver_setup!(plant::PlantState)
    N   = length(plant.statevec)
    lb  = zeros(N)
    ub  = fill(1e9, N)
    A   = sparse(Diagonal(ones(N)))
    qp  = quadratic_problem(plant)
    OSQP.setup!(plant.solver, P=qp.P, q=qp.q, l=lb, u=ub, A=A, verbose=false)
    return qp.P
end

function quadratic_problem(plant::PlantState)
    measurements  = plant.measurements
    
    R⁻¹ = Diagonal(inv.(variance_vector(measurements)))
    P⁻¹ = sparse(plant.stateinv)
    x   = copy(plant.statevec)
    y   = measurement_vector(measurements)

    yh = prediction_vector(measurements, x)
    C  = observation_matrix(measurements, x)

    ε  = (yh - C*x)
    yc = y - ε

    return (
        P = (C'*R⁻¹*C + P⁻¹),
        q = -(C'*R⁻¹*yc + P⁻¹*x)
    )
end

valuelength(x::AbstractVector{<:AbstractMeas}) = sum(valuelength, x, init=0)
valuelength(c::MeasCollection) = sum(fn-> valuelength(c[fn]), fieldnames(MeasCollection))

function observation_matrix(c::MeasCollection, x::AbstractVector{T}) where T
    obsmat = zeros(promote_type(T,Float64), valuelength(c), length(x))
    
    #Fine-tuned shortcut to be implemented later in replacement of jacobian!
    #ForwardDiff.jacobian!(obsmat, z->prediction_vector(c, z), x)
    _fill_array!(addobserver!, obsmat, c, x)

    return sparse(obsmat)
end

variance_vector(c::MeasCollection) = variance_vector!(zeros(Float64, valuelength(c)), c)
variance_vector!(measvec::AbstractVector, c::MeasCollection) = _fill_array!(setvariance!, measvec, c)

measurement_vector(c::MeasCollection) = measurement_vector!(zeros(Float64, valuelength(c)), c)
measurement_vector!(measvec::AbstractVector, c::MeasCollection) = _fill_array!(setmeasurement!, measvec, c)

prediction_vector(c::MeasCollection, x::AbstractVector{T}) where T = prediction_vector!(zeros(T, valuelength(c)), c, x)
prediction_vector!(measvec::AbstractVector, c::MeasCollection, x::AbstractVector) = _fill_array!(setprediction!, measvec, c, x)


function _fill_array!(fillfunc!::Function, arr::AbstractArray, c::MeasCollection, args...)
    function multi_fill!(arr, ind, vmeas, args...)
        for meas in vmeas
            fillfunc!(arr, ind, meas, args...)
        end
        return nothing
    end

    #Mapping on a tuple removes type-instability
    ind = Ref(0)
    map(fn-> multi_fill!(arr, ind, c[fn], args...), fieldnames(MeasCollection))

    return arr
end

function state_indices(m::MoleBalance)
    indices = Int64[]
    for inlet in m.inlets
        state_indices!(indices, inlet)
    end
    for outlet in m.outlets
        state_indices!(indices, outlet)
    end
    return sort!(indices)
end
state_indices(m::AbstractMeas) = state_indices(m.stream)


function state_indices!(indices::AbstractVector, stream::StreamRef)
    union!(indices, stream.comp)
    if hasparents(stream)
        push!(indices, stream.flow)
    end
    return indices
end
state_indices(stream::StreamRef) = sort!(state_indices!(Int64[], stream))


function addobserver!(obsmat::AbstractMatrix, rows::RefValue{<:Integer}, m::AbstractMeas,  statevec::AbstractVector)
    obsview = ObserverView(m::AbstractMeas, statevec)
    addobserver!(obsmat, nextindex!(rows, m), obsview)
    return obsmat 
end

function addobserver!(obsmat, rows::Integer, obsview::ObserverView)
    obsmat[rows:rows, obsview.stateind] += obsview.jacobian
end

function addobserver!(obsmat, rows::AbstractVector{<:Integer}, obsview::ObserverView)
    obsmat[rows, obsview.stateind] += obsview.jacobian
end



#=
function reconcile_statevec!(plant::PlantState)
    negloglik(x) = -loglik(x)
    objfunc = OptimizationFunction(negloglik, AutoZygote())
    N = length(plant.statevec)
    lb = fill(1e-6, N)
    ub = fill(1e12, N)

    problem = Optimization.OptimizationProblem(objfunc, copy(plant.statevec), plant, lb=lb, ub=ub)
    results = solve(problem, Optimization.LBFGS(), reltol=1e-6)
    u = clamp.(results.u, lb, ub)

    if loglik(u, plant) > loglik(plant.statevec, plant)
        plant.statevec .= results.u
    end
    return results
end
=#

#=
function observation_matrix(meas::AbstractMeas, xref::AbstractVector)
    obsfunc(x::AbstractVector) = -innovation(x, meas)
    return Zygote.jacobian(obsfunc, xref)[1]
end
=#

function update_balance_errors!(plant::PlantState{S}) where S
    molebalances = plant.measurements.MoleBalance

    for ii in eachindex(molebalances)
        balance = molebalances[ii]
        balance.value = Species{S}(balance.value .- prediction(plant.statevec, balance))
    end

    return plant 
end




