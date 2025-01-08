include("_PlantState.jl")
#using Optim, LineSearches
using Optimization
import Zygote

function errorgradient(statevec::AbstractVector, plant::PlantState)
    return Zygote.gradient(x->negloglik(x, plant), statevec)[1]
end

function reconcile!(plant::PlantState)
    optimresults = reconcile_statevec!(plant)
    statecov = reconcile_statecov!(plant)
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




