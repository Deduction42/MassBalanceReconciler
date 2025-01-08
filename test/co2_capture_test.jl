using Revise
using MassBalanceReconciler
using JSON3
using Dates

#===========================================================================================================
ToDo:

(1) Ensure that the dictionary has a "_UNIX_TIMESTAMP" entry
(2) Remove the "timestamp" argument from readvalues() and use "_UNIX_TIMESTAMP" instead
(3)  -  Ensure that the evaluation interval is properly calculated and updated for all balance measurements
(4) Wrap entire procedure into a reconcile!(plantstate, measdict) call
     -  readvalues!
     -  updatethermo!
     -  predict!
     -  reconcile_statevec!
     -  reconcile_statecov!
     -  update_balance_errors!

===========================================================================================================#


capture_dict = open(joinpath(@__DIR__,"co2_capture_plant.json")) do fh
    JSON3.read(fh, numbertype=Float64)
end

meas_dict = open(joinpath(@__DIR__,"co2_capture_meas.json")) do fh
    meas_json = JSON3.read(fh)
    Dict{String,Float64}(String(k)=>v for (k,v) in pairs(meas_json))
end

plantinfo = PlantInfo(capture_dict)
plantstate = PlantState(plantinfo)

t0 = datetime2unix(floor(now(), Hour(1)))
t1 = t0 + plantstate.clock.interval[]
plantstate.clock.timestamp[] = t0
meas_dict[MassBalanceReconciler.TIMESTAMP_KEY] = t1

#=
readvalues!(plantstate, meas_dict)
updatethermo!(plantstate)
predict!(plantstate, 60)

P0 = deepcopy(plantstate.statecov)
@time optimresults = reconcile_statevec!(plantstate)
P1 = reconcile_statecov!(plantstate)
update_balance_errors!(plantstate)
=#
(optimresults, statecov) = reconcile!(plantstate, meas_dict)

