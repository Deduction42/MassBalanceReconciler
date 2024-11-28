include("_AbstractMeas.jl")
using LinearAlgebra

@kwdef struct PlantState{L, N}
    timestamp    :: Base.RefValue{Float64}
    statevec     :: Vector{Float64}
    statecov     :: Matrix{Float64}
    dpredictor   :: @NamedTuple{A::Matrix{Float64}, Q::Matrix{Float64}}
    measurements :: MeasCollection{L, Float64, N}
    streams      :: Vector{StreamInfo{L,N}}
    nodes        :: Vector{NodeInfo{L,N}}
end

function PlantState(plant::PlantInfo{Lc,Nc}, thermo::Dict{Symbol, <:ThermoState{Ls,<:Real,Ns}}) where {Lc,Nc,Ls,Ns}
    #Read the plant info and re-index it, finding the total state length
    Nx = stateindex!(plant)
    statevec = zeros(Float64, Nx)

    #Get the default flow rate of all streams
    stream_defaults = Dict{Symbol, Species{Lc,Float64,Nc}}()

    for stream in plant.streams
        state = thermo[stream.id]

        #Calculate the molar flow of species
        xs = state.n[:]./sum(state.n[:])
        MW = dot(molar_weights(state)[:], xs)
        nf = stream.massflow/MW
        ns = Species{Ls, Float64, Ns}(nf.*xs)

        #Aggregate the species as components (which may be different from thermodynamic model)
        streamval = total(Species{Lc}, ns)

        #Store results
        stream_defaults[stream.id] = streamval
        statevec[stream.index] = streamval
    end

    #Get the reaction defaults based off limiting reagents of the inputs and assign it to the state
    for node in plant.nodes
        if !isempty(node.reactions)
            inputs = sum(Base.Fix1(speciesvec, statevec), node.inlets)

            for rxn in node.reactions
                extent = stoich_extent(rxn.stoich, inputs)
                statevec[rxn.extent] = extent
            end
        end
    end

    #Build the state covariance assuming the nominal values are the standard deviation
    avgstatevec = deepcopy(statevec)
    for stream in plant.streams
        ind = stream.index
        avgstatevec[ind] .= sum(avgstatevec[ind])/Nc
    end
    statecov = Matrix(Diagonal(avgstatevec.^2))

    #Build the predictor based on relationships, and the noise intensity based off initial state covariance
    A = zeros(Nx,Nx)
    stream_dict = Dict{Symbol, StreamInfo{Lc,Nc}}(x.id=>x for x in plant.streams)
    for relationship in plant.relationships
        streamind = stream_dict[relationship.id].index.data 
        parentind = stream_dict[relationship.parent].index.data 
        
        for (s,p) in zip(streamind, parentind)
            A[s,s] = -1/relationship.timeconst
            A[s,p] = relationship.factor/relationship.timeconst
        end
    end

    #Noise intensity is assumed to be quite large (trust measurements) such that it reaches typical magnitude in 60 seconds
    Q = statecov./60 
    dpredictor = (A=A, Q=Q)


    #Build the measurements based off the thermodynamic information
    meascollection = MeasCollection{Lc,Float64}()

    for measinfo in plant.measurements
        type = measinfo.type
        meas = build(type, measinfo, stream_dict, thermo)
        push!(meascollection[type], meas)
    end

    for nodeinfo in plant.nodes
        meas = build(MoleBalance, nodeinfo, stream_dict)
        push!(meascollection[MoleBalance], meas)
    end


    #Populate the final object with constructed values and pass through the stream and node information
    return PlantState{Lc,Nc}(
        timestamp = Ref(0.0),
        statevec = statevec,
        statecov = statecov,
        dpredictor = dpredictor,
        measurements = meascollection,
        streams = plant.streams,
        nodes = plant.nodes
    )
end

function readvalues!(plant::PlantState, d::AbstractDict, t::Real)
    interval = max(0.0, t - plant.timestamp[])
    readvalues!(plant.measurements, d)
    setintervals!(plant.measurements.MoleBalance, interval)
    return plant
end

function updatethermo!(plant::PlantState, d::Dict{Symbol, <:ThermoState})
    updatethermo!(plant.measurements, d)
    return plant 
end

#=============================================================================
Populate the tag dictionary with translated anlyzer values
=============================================================================#
function translate!(tagdict::Dict{String}, measurements::MeasCollection{L}, thermo::Dict{Symbol, <:ThermoState}) where {L}
    for meas in measurements.MoleAnalyzer
        molefracs  = total(Species{L}, thermo[meas.streamid].n)[:]
        molefracs  = molefracs./sum(molefracs)

        for ii in eachindex(molefracs)
            tagdict[meas.tag[ii]] = molefracs[ii]
        end
    end
end

function predict!(plant::PlantState, interval::Real)
    Ad = exp(interval.*plant.dpredictor.A)
    plant.statevec .= Ad*plant.statevec
    plant.statecov .= Ad'*plant.statecov*Ad + interval.*plant.dpredictor.Q
    return plant 
end

function negloglik(x::AbstractVector, plant::PlantState)
    Δx = x .- plant.statevec
    state_negloglik = Δx'*plant.statecov*Δx
    meas_negloglik  = negloglik(x, plant.measurements)
    return state_negloglik + meas_negloglik
end





