include("_PlantState.jl")
#using Optim, LineSearches
using Optimization
import Zygote

function errorgradient(statevec::AbstractVector, plant::PlantState)
    return Zygote.gradient(x->negloglik(x, plant), statevec)[1]
end


function reconcile!(plant::PlantState, data::AbstractDict{K, TimeSeries{T}}) where {K<:AbstractString, T<:Real}
    #Average the data over the intervals
    data_avg = si_time_averages(plant, data)
    vt = timestamps(first(data_avg))

    #Initialize the samples and states
    samples = Dict{String,Float64}()
    states = TimeSeries([TimeRecord(t0, deepcopy(plant.statevec))])
    stdevs = TimeSeries([TimeRecord(t0, sqrt.(diag(plant.statecov)))])

    for t in vt
        interpolate!(samples, data_avg, t)
        reconcile!(plant, samples)
        push!(states, TimeRecord(t, deepcopy(plant.statevec)))
        push!(stdevs, TimeRecord(t, sqrt.(diag(plant.statecov))))
    end

    return PlantSeries(plant, states, stdevs)
end



function si_time_averages(plant::PlantState, data::AbstractDict{<:String, TimeSeries{T}}) where T
    #Create a vector of sampled timestamps
    Δt = Second(round(plant.clock.interval[]))
    t0 = floor(minimum(x->datetime(x[begin]), values(data_avg)), Δt)
    tN = ceil(maximum(x->datetime(x[end]), values(data_avg)), Δt)
    vt = t0:Δt:tN

    #Obtain tags and initialize the averages
    tags = gettags(plant.measurements)
    data_avg = Dict{String, TimeSeries{promote_type(T,Float64)}}()

    for k in tags #Calculate the averages and convert to SI units
        vec_avg  = average(data[k], vt)
        vec_avg .= to_si_units.(values(vec_avg), plant.units[k])
        data_avg[k] = vec_avg
    end

    return data_avg
end

to_si_units(x::Real, u::Quantity)   = ustrip(uexpand(x*u))
from_si_units(x::Real, u::Quantity) = ustrip(uconvert(u, x*Quantity(1, dimension(uexpand(u)))))

function interpolate!(samples::AbstractDict{K, <:Real}, data::AbstractDict{K, TimeSeries{T}}, t::Real) where {K,T<:Real}
    f = TimeRecords._interpolate_linsat
    indhint = max(findbounds(ts, t)[1], 1)
    for k in keys(data)
        samples[k] = interpolate(f, data[k], t, indhint)
    end
    samples[TIMESTAMP_KEY] = t 
    return samples
end


function reconcile!(plant::PlantState, data::AbstractDict{<:String, <:Real})
    readvalues!(plant, data)
    updatethermo!(plant)
    predict!(plant)
    optimresults = reconcile_statevec!(plant)
    statecov = reconcile_statecov!(plant)
    update_balance_errors!(plant)
    return (optimresults, statecov)
end

function reconcile_statevec!(plant::PlantState)
    objfunc = OptimizationFunction(negloglik, AutoZygote())
    N = length(plant.statevec)
    lb = fill(1e-12, N)
    ub = fill(Inf, N)
    problem = Optimization.OptimizationProblem(objfunc, plant.statevec.*1, plant, lb=lb, ub=ub)
    results = solve(problem, Optimization.LBFGS(), reltol=1e-9)

    if results.objective < negloglik(plant.statevec, plant)
        plant.statevec .= results.u
    end
    return results
end



function reconcile_statecov!(plant::PlantState)
    for measvec in plant.measurements[:]
        for meas in measvec
            reconcile_statecov!(plant, meas)
        end
    end
    return plant.statecov
end


function reconcile_statecov!(plant, meas)
    inds = falses(length(plant.statevec))
    addinds!(inds, meas)

    C  = observation_matrix(meas, plant.statevec)[:,inds]
    R  = noisecov(meas)
    Pi = plant.statecov[inds, inds]
    S  = C*Pi*C' .+ R
    K  = (Pi*C')/S
    Kp = I-K*C

    if all(isfinite, K)
        #Update inner covariance structure
        plant.statecov[inds, inds] .= Kp*Pi*Kp' .+ K*R*K'

        #Get the outer indices
        outs = .!inds

        #Update the outer covariance
        PU = Kp*(@view plant.statecov[inds, outs])
        plant.statecov[inds, outs] .= PU
        plant.statecov[outs, inds] .= PU'
    end

    return plant.statecov
end

function observation_matrix(meas::AbstractMeas, xref::AbstractVector)
    obsfunc(x::AbstractVector) = -innovation(x, meas)
    return Zygote.jacobian(obsfunc, xref)[1]
end

function update_balance_errors!(plant::PlantState{S}) where S
    molebalances = plant.measurements.MoleBalance

    for ii in eachindex(molebalances)
        balance = molebalances[ii]
        updated = (@set balance.value = Species{S}(balance.value .- prediction(plant.statevec, balance)))
        molebalances[ii] = updated
    end

    return plant 
end




