using JSON
include(joinpath(@__DIR__, "Assembly.jl"))

@time assetList = JSON.parsefile("TestSystem.json")
@time measList  = JSON.parsefile("TestSystemMeas.json")

@time assetBOM   = build_bom(assetList)
@time measGroups = collect_by_type( values(build_bom(measList)) )

@time bomBounds = build_bom(assetList, elType=Bounds)
@time (indexer, n) = bom_indexer(assetBOM)

bndVec = get_bom_state!(fill(Bounds(NaN,NaN),n), bomBounds, indexer)

LB = getproperty.(bndVec, :LB)
UB = getproperty.(bndVec, :UB)

res = zeros( sum(length.(values(measGroups))) )
@time predict_measurements!(res, measGroups, assetBOM)
