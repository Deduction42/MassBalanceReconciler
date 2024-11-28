include("_PlantState.jl")
using Optim, LineSearches
import Zygote

function errorgradient(statevec::AbstractVector, plant::PlantState)
    return Zygote.gradient(x->negloglik(x, plant), statevec)[1]
end

function reconcile!(plant::PlantState)
    optimresults = reconcile_statevec!(plant)
    return optimresults
end


function reconcile_statevec!(plant::PlantState)
    objfunc(x::AbstractVector) = negloglik(exp.(x), plant)

    function objgrad!(g::AbstractVector, x::AbstractVector) 
        g .= Zygote.gradient(objfunc, x)[1]
        return g
    end

    initial = log.(max.(1e-3, plant.statevec))

    # requires using LineSearches
    inner_optimizer = LBFGS(linesearch=LineSearches.HagerZhang())
    options = Optim.Options(f_reltol=1e-9, g_tol=0)
    results = optimize(objfunc, objgrad!, initial, inner_optimizer, options)

    #Overwrite the state if successful
    if objfunc(results.minimizer) < objfunc(log.(plant.statevec))
        plant.statevec .= exp.(results.minimizer)
    end

    return results
end

function reconcile_statecov!(plant)
    for measvec in plant.measurements[:]
        for meas in measvec
            reconcile_statecov!(plant, meas)
        end
    end
    return plant.statecov
end


function reconcile_statecov!(plant, meas)
    inds = stateindex(meas)
    C  = observation_matrix(meas, plant.statevec)[:,inds]
    R  = noisecov(meas)
    Pi = plant.statecov[inds, inds]
    S  = C*Pi*C' .+ R
    K  = (Pi*C')/S
    Kp = I-K*C

    display(meas.id)
    display(Kp)


    if all(isfinite, K)
        #Update inner covariance structure
        plant.statecov[inds, inds] .= Kp*Pi*Kp' .+ K*R*K'

        #Get the outer indices
        outbool = trues(length(plant.statevec))
        outbool[inds] .= false
        outs = findall(outbool)

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

#=
#Older method that uses box constraints
function reconcile!(plant::PlantState)
    objfunc(x::AbstractVector) = negloglik(x, plant)

    function objgrad!(g::AbstractVector, x::AbstractVector) 
        g .= Zygote.gradient(objfunc, x)[1]
        return g
    end

    lower = fill(zero(Float64), length(plant.statevec))
    upper = fill(Inf64, length(plant.statevec))
    initial = max.(1e-9, plant.statevec)

    # requires using LineSearches
    inner_optimizer = LBFGS(linesearch=LineSearches.HagerZhang())
    options = Optim.Options(f_reltol=1e-6)
    results = optimize(objfunc, objgrad!, lower, upper, initial, Fminbox(inner_optimizer), options)

    #Overwrite the state if successful
    if objfunc(results.minimizer) < objfunc(plant.statevec)
        plant.statevec .= results.minimizer
    end

    return results
end
=#