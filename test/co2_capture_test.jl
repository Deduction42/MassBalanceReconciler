using Revise
using MassBalanceReconciler
using JSON3

capture_dict = open(joinpath(@__DIR__,"co2_capture_plant.json")) do fh
    JSON3.read(fh, numbertype=Float64)
end

meas_dict = open(joinpath(@__DIR__,"co2_capture_meas.json")) do fh
    meas_json = JSON3.read(fh)
    Dict{String,Float64}(String(k)=>v for (k,v) in pairs(meas_json))
end

plantinfo = PlantInfo(capture_dict)
plantstate = PlantState(plantinfo)

readvalues!(plantstate, meas_dict, 60*15)
updatethermo!(plantstate)
predict!(plantstate, 60)

P0 = deepcopy(plantstate.statecov)
@time optimresults = reconcile_statevec!(plantstate)
P1 = reconcile_statecov!(plantstate)
update_balance_errors!(plantstate)


