###############################
# Epistasis functions #########
# #############################

function computeDoubleMutGSE(Q::Network,
                             singleMutStructs::Matrix{Matrix{Matrix{Float64}}},
                             ligDex::Vector{CartesianIndex{2}};
                             numTrials::Int=1000)

    # compute the ground state energies for all
    # double mutants. Stores them in an Vector of Matricies.

    function kernel(Q::Network, indices::Tuple{Int, Int},
                    structs::Vector{Vector{Matrix{Matrix{Float64}}}},
                    ligDex, numTrials)
        # warning This function only work when numtypes = 5
        
        Q = deepcopy(Q)
        i,j = indices
        structs_i,  structs_j = structs
        numTypes = 5
        numMuts = numTypes - 1 # number of mutations to make
        numLigs = length(Q.ligs)
        out = Array{Array{Float64,2},2}(undef, numMuts, numMuts)
        f(i) = setdiff( 1:numTypes, Q.geno[i]) 
        availTypes_i = f(i)
        availTypes_j = f(j)
        energies = zeros(numLigs, numLigs)
        for m in 1:numMuts, n in 1:numMuts
            ti = availTypes_i[m]
            tj = availTypes_j[n]
            mutateNodeCouple!(Q, i, j, ti, tj)
            reg = [structs_i[ti][ligDex]; structs_j[tj][ligDex] ] 
            energies[ligDex], _ = findGroundStates(Q, reg, ligDex, numTrials)
            out[m,n] = copy(energies)
        end
        return out
    end

    Q = deepcopy(Q) # to detach from Q outside function
    numNodes = size(Q.A, 1)
    numTypes = size(Q.typeTable,1)

    GSE_WT = Q.energies
    # only mutate up to five types
    if numNodes > 5; num2Mut=5; else; num2Mut=numNodes; end

    # build input
    N = Int( 0.5*(numNodes^2 - numNodes)) # num element in lower tri
    indices = Vector{Tuple{Int64,Int64}}(undef, N)
    structsList = Vector{Vector{Vector{Matrix{Matrix{Float64}}}}}(undef,N)
    k = 0
    for j in 1:numNodes-1, i in j+1:numNodes
        k += 1
        indices[k] = (i,j)
        structsList[k] = [ singleMutStructs[i,:], singleMutStructs[j,:] ]
    end

    # compute double mutant GSE
    return pmap( (x,y) -> kernel(Q, x, y, ligDex, numTrials), indices, structsList)
end



function doubleMutGSE2Fit(dmGSE::Vector{Matrix{Matrix{T}}},
                          params::Params) where T<:Number
    # convert double mutant GSE array to double mutant
    # fitness array
    f(m) = map(x -> computeFitness(x, params), m)
    return f.(dmGSE)
end


function computeEpistasis(f0::Real,  # WT fitness
                          f1::Real,  # single mutant 1 fitness
                          f2::Real,  # single mutant 2 fitness 
                          f12::Real) # double mutant fitness
    return f12 - f1 - f2 + f0
end


function computeEpistasisArray(geno_WT::Vector{Int},
                               fit_WT::T,
                               fit_singleMuts::Matrix{T},
                               fit_doubleMuts::Vector{Matrix{T}}
                              ) where T<:Number
    
    # compute Epistatis terms. These terms only apply to double mutants.
    # Return output in linear lower diagonal form.

    Epi = deepcopy(fit_doubleMuts)
    numNodes, numTypes  = size(fit_singleMuts)
    tmp = zeros(numTypes-1, numTypes-1)

    for k in eachindex(fit_doubleMuts)
        i,j = k2ij(k,numNodes) 
        typei = geno_WT[i]
        typej = geno_WT[j]
        for (m,typem) in enumerate(setdiff(1:numTypes,typei))
            fit_sm_1 = fit_singleMuts[i,typem] 
            for (n,typen) in enumerate(setdiff(1:numTypes,typej))
                fit_sm_2 = fit_singleMuts[j,typen]
                fit_dm = fit_doubleMuts[k][m,n]
                Epi[k][m,n] = computeEpistasis(fit_WT, 
                                               fit_sm_1,
                                               fit_sm_2,
                                               fit_dm)
            if i==17 && j==1
            #    println("$typei $typej")
            #    println(fit_WT)
            #    println(fit_sm_1)
            #    println(fit_sm_2)
            #    println(fit_dm)
            #    println(Epi[k][m,n])
            end


            end
        end
    end
    return Epi
end


function makeEpiMat(df::Matrix{T},
                    Epi::Vector{Matrix{T}}
                   ) where T<:Number
    # df is the matrix of mutational sensitivities
    # Epi is the Epistatic array
    
    numNodes, numTypes = size(df)
    epiMat = zeros(numNodes,numNodes)
    f(mat) = mean(abs.(mat)) # this might need to be changed.
   # f(mat) = mean(mat) # this might need to be changed.
    vec2UpperLowerTriMat!(epiMat, f.(Epi))

    epiMat[diagind(epiMat)] .= vec(mean(abs.(df), dims=2)) .* (numTypes/ (numTypes-1))
    return epiMat
end



function fastTrackEpiMat(Q::Network,
                         singMutGSE::Matrix{Matrix{Float64}},
                         doubMutGSE::Vector{Matrix{Matrix{Float64}}},
                         assay::String)

    
    # convert ground state energies to fitnesses
    params = constructParams(Q, assay)
    fit_WT = computeFitness(Q, assay)
    singMutFit = map(x -> computeFitness(x, params), singMutGSE)
    doubMutFit = doubleMutGSE2Fit(doubMutGSE, params)
    
    # convert fitnesses to epistatic terms
    diagEpiTerms = singMutFit .- fit_WT
    offDiagEpiTerms = computeEpistasisArray(Q.geno, fit_WT, singMutFit, doubMutFit)

    epiMat = makeEpiMat(diagEpiTerms, offDiagEpiTerms)
    return epiMat, offDiagEpiTerms
end



#function ensembleComputeDoubleMutGSE(ensemble::Array{Network},
#                                     singleMutStructs::Matrix{Matrix{Matrix{Matrix{Float64}}}}, # ens, dms, 3x3, xy
#                                     ligDex::Vector{CartesianIndex{2}};
#                                     numTrials::Int=1000)
#
#    
#
#end









