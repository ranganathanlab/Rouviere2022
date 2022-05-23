function computeGSEMutSens(Q::Network,
                           ligDex::Vector{CartesianIndex{2}},
                           mutType::String;
                           GS_settings::Tuple{Int, Int}=(10,100),
                           saveStructs::Bool=false,
                           energyThresh::Number=1e-5,
                           bondStrain::Number=0.2)

    possibleMutTypes = ["Node", "Bond", "SurfaceBond", "NodeCouple"]
    @assert mutType in possibleMutTypes

    Q = deepcopy(Q) # to detach from Q outside function
    Qwildtype = deepcopy(Q)
    GSE_WT = copy(Q.energies)

    
    # Get information on how to do the mutational scan for the type of mutant pertubation.
    mutDetailsList = getMutDetailsList(Q, mutType; bondStrain)
    # Do the deep mutational scan
    GSE_DMS, structs_DMS, registryLong =  doMutScan(Q, mutType, mutDetailsList, ligDex, GS_settings)

    # remove redundancy in the list of saved structure
    removeAllLig!(Q)
    registry = removeRedundantStructures(Q, registryLong)

    # test wildtype with new states
    GSE_WT[ligDex], structs_WT = findGroundStates(Q, ligDex, registry; energyThresh)

    # test all mutants with new states
    testAllMutantsWithAllStructures!(GSE_DMS, structs_DMS, Q, mutType, mutDetailsList, ligDex, registry, energyThresh)

    if saveStructs
        return GSE_WT, GSE_DMS, structs_WT, structs_DMS
    else
        return GSE_WT, GSE_DMS
    end
end

function doMutScan(Q::Network, mutType::String, mutDetailsList, ligDex, GS_settings)
    # HELPER FUNCTION FOR computeGSEMutSens.
    # Do the deep mutational scan saving ground state energies
    # and structures for all mutants.

    # intialize data structures to store ground state energy and structures 
    Q = deepcopy(Q) # detach from outside
    Qwildtype = deepcopy(Q)

    numLigs = length(Q.ligs)
    numNodes = size(Q.A, 1)
    GSE_DMS = Matrix{Matrix{Float64}}(undef, size(mutDetailsList))
    [GSE_DMS[i] = copy(Q.energies) for i in eachindex(GSE_DMS)] 
    structs_DMS = Matrix{Matrix{Matrix{Float64}}}(undef, size(mutDetailsList))
    [structs_DMS[i] = deepcopy(Q.structs) for i in eachindex(structs_DMS)]

   # structs_packet = [zeros(numNodes,2) for i in 1:numLigs, j in 1:numLigs]
   # [structs_DMS[i] = deepcopy(structs_packet) for i in eachindex(structs_DMS)]
    registryLong = Matrix{Float64}[]


    for (i, mutDetails) in enumerate(mutDetailsList)
        mutate!(Q, mutType, mutDetails)
        GSE_DMS[i][ligDex], structs_DMS[i][ligDex] = findGroundStates(Q, ligDex, GS_settings)
        
        append!(registryLong, structs_DMS[i][ligDex]) 
        # undo the mutation
        Q = deepcopy(Qwildtype)
    end
    return GSE_DMS, structs_DMS, registryLong
end

function testAllMutantsWithAllStructures!(GSE_DMS,
                                          structs_DMS,
                                              Q::Network,
                                          mutType::String,
                                          mutDetailsList,
                                          ligDex,
                                          registry,
                                          energyThresh)
    # test all mutants with all of the unique ground state structures.
    Q = deepcopy(Q)
    Qwildtype = deepcopy(Q)
    for (i, mutDetails) in enumerate(mutDetailsList)
        mutate!(Q, mutType, mutDetails)
        Q.energies = GSE_DMS[i]
        Q.structs = structs_DMS[i]
        GSE_DMS[i][ligDex], structs_DMS[i][ligDex] = findGroundStates(Q, ligDex, registry; energyThresh)
        Q = deepcopy(Qwildtype)
    end
    return nothing
end

function removeRedundantStructures(Q::Network,
                                   structureList::Vector{Matrix{Float64}};
                                   energyThresh::Number=1e-5)

    refinedStructList = Matrix{Float64}[]
    energyList = Float64[]
    for i in eachindex(structureList)
        e, xyr = relaxSprings(structureList[i], Q.pheno)
        if !intol(e, energyList, energyThresh)
            push!(energyList, e)
            push!(refinedStructList, copy(xyr))
        end
    end
    return refinedStructList
end





function getMutDetailsList(Q::Network,
                           mutType::String;
                           bondStrain::Number=0.2)

    # This function is designed to work with computeGSEMutSens
    # and returns the list of mutational details needed to compute
    # a variety GSE deep mutational scans.

    if mutType == "Node"
        numNodes = length(Q.geno)
        numTypes =  size(Q.typeTable,1)
        
        numPossMuts = numTypes - 1
        numPossMuts > 4 ? num2Mut = 4 : num2Mut = numPossMuts # only make 4 mutations max.
        mutDetailsList = Matrix{Any}(undef, numNodes, num2Mut)
        for i in 1:numNodes
            listPossMutants = setdiff(1:numTypes, Q.geno[i])
            for j in 1:num2Mut
                mutDetailsList[i,j] = (i,listPossMutants[j])
            end
        end

    elseif mutType == "Bond"
        error("Bond option is not set up yet!!")
    elseif mutType == "SurfaceBond"
        bondStrain = abs(bondStrain)
        bondStrains = [bondStrain, -bondStrain]
        surfaceBonds = getSurfaceBonds()
        mutDetailsList = Matrix{Any}(undef, length(surfaceBonds), length(bondStrains))
        for i in 1:length(surfaceBonds), j in 1:length(bondStrains)
            mutDetailsList[i,j] = (surfaceBonds[i], bondStrains[j])
        end

    elseif mutType == "NodeCouple"
        error("NodeCouple option is not set up yet!!")
    else
        error("Your mutType was not one of the following: Node, Bond, SurfaceBond, NodeCouple")
    end
    return mutDetailsList
end

function getSurfaceBonds()
        # these are correct only for my 49 node network.
    return [CartesianIndex(22,29),
            CartesianIndex(29,36),
            CartesianIndex(36,43),
            CartesianIndex(43,44),
            CartesianIndex(44,1),
    #        CartesianIndex(1,7),
            CartesianIndex(7,48),
            CartesianIndex(48,49),
            CartesianIndex(49,42),
            CartesianIndex(42,35),
            CartesianIndex(35,28)]
end

    
function computeSingleMutsFit(GSE_DMS::Array{<:Array{Float64}}, # either Mut, Bond, Allo
                              params::Params)
    # compute the fitness of each single mutant in the DMS
    return map(x -> computeFitness(x, params), GSE_DMS )
end

function computeFitSens(GSE_WT::T,
                        GSE_DMS::Array{T}, # either Mut, Bond, Allo
                        params::Params) where T <: Array{Float64}
    # compute the mutational effect of each mutant in the DMS
    return computeSingleMutsFit(GSE_DMS, params) .- computeFitness(GSE_WT, params)
end


function computeSurfaceAllostery(GSE_WT::Matrix{Float64},
                          GSE_DMS::Array{<:Matrix{<:Float64}},
                          params::Params,
                          coopSign::String)::Float64
    # compute a measure of allostery 
    df = computeFitSens(GSE_WT, GSE_DMS, params)
    if coopSign == "+"
        return maximum(df)
    elseif coopSign == "-"
        return minimum(df)
    else
        error("sign must be \"+\" or \"-\"")
    end
    return nothing
end

function findAllostericHotSpot(df::Matrix{Float64}, coopSign::String)
    if coopSign == "+"
        index = argmax(df)
    elseif coopSign == "-"
        index = argmin(df)
    else
        error("sign must be \"+\" or \"-\"")
    end
    surfaceBonds = getSurfaceBonds()
    hotSpot = surfaceBonds[index[1]]
    bondStrainSign = -2*index[2] +3
    return hotSpot, bondStrainSign
end
































