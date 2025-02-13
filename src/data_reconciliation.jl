include("_PlantState.jl")
#using Optim, LineSearches
using Optimization
import Zygote

function errorgradient(statevec::AbstractVector, plant::PlantState)
    return Zygote.gradient(x->negloglik(x, plant), statevec)[1]
end


function reconcile!(plant::PlantState, interval::TimeInterval, data::AbstractDict{K, TimeSeries{T}}) where {K<:AbstractString, T<:Real}
    #Average the data over the intervals
    data_avg = si_time_averages(plant, interval, data)
    vt = data_avg[TIMESTAMP_KEY]

    #Initialize the samples and states
    samples = Dict{String,Float64}()
    t0 = interval[begin]
    states  = TimeSeries([TimeRecord(t0, deepcopy(plant.statevec))])
    stdevs  = TimeSeries([TimeRecord(t0, sqrt.(diag(plant.statecov)))])

    #Reconciliation for each timestamp
    for ii in eachindex(vt)
        for (k, v) in pairs(data_avg)
            samples[k] = v[ii]
        end
        t  = samples[TIMESTAMP_KEY]
        optimresults = reconcile!(plant, samples)[1]

        @info "$(round(ii/length(vt)*100)) % complete 
         status: $(optimresults.retcode)
         time: $(round(optimresults.stats.time, digits=6)) s
         iterations: $(optimresults.stats.iterations)"

        push!(states, TimeRecord(t, deepcopy(plant.statevec)))
        push!(stdevs, TimeRecord(t, sqrt.(diag(plant.statecov))))
    end

    return PlantSeries(plant, states, stdevs)
end


"""
si_time_averages(plant::PlantState, interval::TimeInterval, data::AbstractDict{<:String, TimeSeries{T}}) where T

Use numerical integration to find averages on 'data', over the time interval 'interval' at the desired sampling rate 
defined in 'plant.clock.interval'. This produces a regularly sampled Dict{String, Vector{T}} where values are in SI units
"""
function si_time_averages(plant::PlantState, interval::TimeInterval, data::AbstractDict{<:String, TimeSeries{T}}) where T
    #Create a vector of sampled timestamps
    Δt = round(plant.clock.interval[])
    vt = interval[begin]:Δt:interval[end]

    #Obtain tags and initialize the averages
    data_avg = Dict{String, Vector{promote_type(T,Float64)}}()

    for taginfo in plant.tags #Calculate the averages and convert to SI units
        vec_avg  = values(average(data[taginfo.tag], vt, order=0))
        vec_avg .= to_si_units.(vec_avg, taginfo.units)
        data_avg[taginfo.tag] = vec_avg
    end

    #Timestamps representing the end of the averaging period
    data_avg[TIMESTAMP_KEY] = vt[(begin+1):end]

    return data_avg
end

to_si_units(x::Real, u::Quantity)   = ustrip(uexpand(x*u))
from_si_units(x::Real, u::Quantity) = ustrip(uconvert(u, x*Quantity(1, dimension(uexpand(u)))))


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
    
    K = try
        (Pi*C')/S
    catch err
        @warn "gain for $(meas.id) produced $(string(err)), skipping covariance update"
        return plant.statecov
    end

    if all(isfinite, K)
        Kp = I-K*C

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




