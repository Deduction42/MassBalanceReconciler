using Revise
using MassBalanceReconciler
using JSON3
using Dates
using TimeRecords

capture_dict = open(joinpath(@__DIR__,"co2_capture_plant.json")) do fh
    JSON3.read(fh, numbertype=Float64)
end

meas_dict = open(joinpath(@__DIR__,"co2_capture_meas.json")) do fh
    meas_json = JSON3.read(fh)
    Dict{String,Float64}(String(k)=>v for (k,v) in pairs(meas_json))
end

#Timestamp definition and construction
t0 = datetime2unix(floor(now(), Hour(1)))
tN = datetime2unix(unix2datetime(t0) + Hour(24))
vt = t0:60:tN
N  = length(vt)
interval = TimeInterval(t0,tN)

#Build timeseries by repeating the same measurements (with noise) from meas_dict
ts_dict = Dict(k=>TimeSeries(TimeRecord.(vt, v.*exp.(0.1.*randn(N)))) for (k,v) in pairs(meas_dict))

#Build plant object, and assign the timestamp
plantinfo  = PlantInfo(capture_dict)
plantstate = PlantState(plantinfo)
plantstate.clock.timestamp[] = t0

#Apply data reconciliation to the timeseries
plantseries = reconcile!(plantstate, interval, ts_dict)