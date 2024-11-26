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
        streamval = Species{Lc}(ns)

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
    statecov = Matrix(Diagonal(statevec.^2))

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

function readvalues!(plant::PlantState, d::AbstractDict)
    readvalues!(plant.measurements, d)
    return plant
end

#=============================================================================
Populate the tag dictionary with translated anlyzer values
=============================================================================#
function translate!(tagdict::Dict{String}, measurements::MeasCollection{L}, thermo::Dict{Symbol, <:ThermoState}) where {L}
    for meas in measurements.MoleAnalyzer
        molefracs  = Species{L}(thermo[meas.streamid].n)[:]
        molefracs  = molefracs./sum(molefracs)

        for ii in eachindex(molefracs)
            tagdict[meas.tag[ii]] = molefracs[ii]
        end
    end
end

#Build the readvalues function for PlantInfo, and ThermoInfo (needed for volume flows)

#Build the update functionality based off the KalmanFilter Dynamic Data Reconciliation model




