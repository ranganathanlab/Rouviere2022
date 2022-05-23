##################################
#### Ground State Functions #######
###################################


function findGroundState(x::Vector{Float64},
                         SM::Matrix{Float64},
                         RLM::Matrix{Float64},
                         EM::Matrix{Float64},
                         settings::Tuple{Int, Int})

    popSize = settings[1]
    numTrials = settings[2]
    p = 6
    X = zeros(length(x), 2popSize)
    E = zeros(2popSize)
    E[1], X[:,1] = relaxSprings(x, SM, RLM, EM)
    x0 = xy2x(buildxy0())
    
    # initial peturb and relax  
    for i in 2:2popSize
        i <= popSize ? x_tmp = x : x_tmp = x0
        xp = applyRandStrain(x_tmp, p)
        E[i], X[:,i] = relaxSprings(xp, SM, RLM, EM)
    end

    # interactively evolve the population
    p = 5
    for j in 1:numTrials
        winners = sortperm(E)[1:popSize]
        losers = sortperm(E)[popSize+1:2popSize]
        for i in 1:popSize
            i_w = winners[i]
            i_l = losers[i]
            xp = applyRandStrain(X[:,i_w],p)
            E[i_l], X[:,i_l] = relaxSprings(xp,SM, RLM, EM)
      end
    end
    e,i = findmin(E)
    return e, X[:,i]
end


function findGroundState(xy::Matrix{Float64},
                         SM::Matrix{Float64},
                         RLM::Matrix{Float64},
                         EM::Matrix{Float64},
                         settings::Tuple{Int, Int})
    E, x = findGroundState(xy2x(xy), SM, RLM, EM, settings)
    return E, x2xy(x)
end



function findGroundStates(Q::Network,
                          ligDex::Vector{CartesianIndex{2}},
                          GS_settings::Tuple{Int, Int})
    
    Q = deepcopy(Q)
    # compute ground state energies for each ligand network complex
    numLigs = length(ligDex)
    energies = zeros(numLigs)
    energies0 = Q.energies[ligDex]
    structs = deepcopy(Q.structs[ligDex])
    structs0 = deepcopy(Q.structs[ligDex])
    pheno0 = deepcopy(Q.pheno)

    # find ground state for each ligand
    for (i,l) in enumerate(ligDex)
        addLig!(Q, l[1], 1) # put ligand at active site
        addLig!(Q, l[2], 2) # put ligand at allosteric site
        energies[i], structs[i] = findGroundState(Q.structs[l], Q.pheno[1], Q.pheno[2], Q.pheno[3], GS_settings)
    end
    Q.pheno = deepcopy(pheno0) # un do binding of ligand


    # Swap ligs and relax
    Q.energies[ligDex] = energies
    Q.structs[ligDex] = structs
    E_mat, confs_mat, trans_mat = swapLigs(Q, ligDex)
    energies = vec(minimum(E_mat, dims=1))
    indices = findLowEIndices(E_mat)
    for i in 1:length(ligDex)
        structs[i] = confs_mat[indices[i],i]
    end

    # set grounds states back to initial.
    Q.structs[ligDex] .= structs0
    Q.energies[ligDex] .= energies0
    return energies, structs
end


function findGroundStates(Q::Network,
                          ligDex::Vector{CartesianIndex{2}},
                          registry::Vector{Matrix{Float64}};
                          energyThresh::Number=1e-5)
    # Find Ground states by relaxing from a resistry of
    # structures, not by the Basin Hopping Alg. 

    energies = copy(Q.energies[ligDex])
    structs = deepcopy(Q.structs[ligDex])
    pheno0 = deepcopy(Q.pheno)
    for (n,l) in enumerate(ligDex)
        addLig!(Q, l[1], 1) # put ligand at active site
        addLig!(Q, l[2], 2) # put ligand at allosteric site
        e = energies[n]
        for k in eachindex(registry)
            e_new, xy_new = relaxSprings(registry[k], Q.pheno)
            if e_new < e - energyThresh
                energies[n] = e =  e_new
                structs[n] = xy_new
            end
        end
    end
    Q.pheno = deepcopy(Q.pheno)
    return energies, structs
end




function findGroundStates!(Q::Network, ligDex::Vector{CartesianIndex{2}},
                           GS_settings::Tuple{Int, Int})
    Q.energies[ligDex], Q.structs[ligDex] = findGroundStates(Q, ligDex, GS_settings)
    return nothing
end



function swapLigs(Q::Network, ligDex::Vector{CartesianIndex{2}})
    # find lowest energy configurations with sol, right, wrong.
    # then swaps ligands and relaxes.
    # returns structures, energies and same state info.

    pheno0 = deepcopy(Q.pheno)
    N = length(ligDex)
    E_mat = zeros(N, N)
    confs_mat = Array{Array{Float64,2},2}(undef,N,N)
    trans_mat = eye(N)
    E_mat[diagind(E_mat)] .= Q.energies[ligDex]
    confs_mat[diagind(confs_mat)] .= Q.structs[ligDex]

    for i in 2:N
        for j in 1:i-1
            # fill lower triangle
            addLig!(Q, ligDex[j][1], 1)
            addLig!(Q, ligDex[j][2], 2)
            E_mat[i,j], confs_mat[i,j] = relaxSprings(Q.structs[ligDex[i]], Q.pheno, gradTol=1e-10)
            trans_mat[i,j] = sameConf(confs_mat[j,j], confs_mat[i,j])
            Q.pheno = deepcopy(Q.pheno)
            
            # fill upper triangle
            addLig!(Q, ligDex[i][1], 1)
            addLig!(Q, ligDex[i][2], 2)
            E_mat[j,i], confs_mat[j,i] = relaxSprings(Q.structs[ligDex[j]], Q.pheno, gradTol=1e-10)
            trans_mat[j,i] = sameConf(confs_mat[i,i], confs_mat[j,i])
            Q.pheno = deepcopy(Q.pheno)
        end
    end
    return E_mat, confs_mat, trans_mat
end





