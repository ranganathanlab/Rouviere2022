
function mutateBond!(Q::Network,
                     bondIndex::CartesianIndex,
                     strain::Float64)
    RLM = Q.pheno[2]
    RLM[bondIndex] = RLM'[bondIndex] *= (1+strain)
    return nothing
end


function mutateNode!(Q::Network,
                     nodeIndex::Int64,
                     newType::Int64)
    numNodes = size(Q.A,1)
    Q.geno[nodeIndex] = newType
    for i in 1:numNodes
        if Q.A[i,nodeIndex] != 0
            typei = Q.geno[i]

            new = Q.typeTable[typei, newType, 1]
            Q.pheno[1][i,nodeIndex] = Q.pheno[1]'[i,nodeIndex] = new

            new = Q.typeTable[typei, newType, 2]
            Q.pheno[2][i,nodeIndex] = Q.pheno[2]'[i,nodeIndex] = new + Q.NRLM[i,nodeIndex] 

            new = Q.typeTable[typei, newType, 3]
            Q.pheno[3][i,nodeIndex] = Q.pheno[3]'[i,nodeIndex] = new

        end
    end
    return nothing
end



function chooseRandomMutation(Q::Network)
    numNodes = size(Q.A, 1)
    numTypes = size(Q.typeTable,1)     
    nodeIndex = rand(1:numNodes)
    newType = rand(setdiff( 1:numTypes, Q.geno[nodeIndex] ))
    return nodeIndex, newType
end


function randomMutateNode!(Q::Network)
    # randomly mutates node
    nodeIndex, newType = chooseRandomMutation(Q)
    mutateNode!(Q, nodeIndex, newType)
    return nothing
end


function mutateNodeCouple!(Q::Network,
                           node1::Int, node2::Int,
                           type1::Int, type2::Int)

    mutateNode!(Q, node1, type1)
    mutateNode!(Q, node2, type2)
    return nothing
end 


function mutate!(Q::Network, mutType::String, mutDetails)
    # Generalized mutatation function that can mutate a 
    # network in many different ways.

    if mutType == "Node"
        mutateNode!(Q, mutDetails...)
    elseif mutType == "Bond" || mutType == "SurfaceBond"
        mutateBond!(Q, mutDetails...)
    elseif mutTupe == "NodeCouple"
        mutateCouple!(Q, mutDetails...)
    end
    return nothing
end








































