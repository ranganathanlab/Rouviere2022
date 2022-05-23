function computeGSEMutSens(Q::Network;
                           ligDex=setLigDex("Big5"),
                           GS_settings::Tuple{Int, Int}=(10,100),
                           saveStructs::Bool=false,
                           energyThresh::Number=1e-5)
    # compute the ground state energies for a network
    # for all single mutants

    Q = deepcopy(Q) # to detach from Q outside function
    numNodes = size(Q.A, 1)
    numTypes = size(Q.typeTable,1)
    numLigs = length(Q.ligs)
    Qwildtype = deepcopy(Q)
    energies = zeros(size(Q.energies))
    GSE_WT = copy(Q.energies)
    slist = Matrix{Float64}[]
   
    # only mutate up to five types
    if numTypes > 5; num2Mut=5; else; num2Mut=numTypes; end
    GSE_DMS = Matrix{Matrix{Float64}}(undef, numNodes, num2Mut)
  
    # initialize data structure to save the network configs 
    structs_DMS = Matrix{Matrix{Matrix{Float64}}}(undef, numNodes, num2Mut)
    structs_packet = [zeros(numNodes,2) for i in 1:numLigs, j in 1:numLigs]
    [structs_DMS[i] = deepcopy(structs_packet) for i in eachindex(structs_DMS)]
   

    # Try to find the ground state energies for all single mutants
    for i in 1:numNodes
        # null mutation
        GSE_DMS[i, Q.geno[i]] = copy(GSE_WT)
        for j in setdiff(1:num2Mut, Q.geno[i])
            mutateNode!(Q, i, j)
            tmpGSE = copy(energies)
            tmpGSE[ligDex], structs = findGroundStates(Q, ligDex, GS_settings)
            GSE_DMS[i,j] = tmpGSE
            structs_DMS[i,j][ligDex] = structs
            append!(slist, structs[:]) 
            # undo the mutation
            Q = deepcopy(Qwildtype)
        end
    end

    # only keep unique structures, compared by their energy. store them in slist_short.
    elist_short = Float64[]
    slist_short = Matrix{Float64}[]
    for i in eachindex(slist)
        e, xyr = relaxSprings(slist[i], Q.pheno)
        if !intol(e, elist_short, energyThresh) # make this a variable
            push!(elist_short, e)
            push!(slist_short, copy(xyr))
        end
    end
    

    # test all single mutants with all of the unique ground state structures.
    for i in 1:numNodes, j in setdiff(1:num2Mut, Q.geno[i])
        mutateNode!(Q, i, j)
        for (n,l) in enumerate(ligDex)
            addLig!(Q, l[1], 1) # put ligand at active site
            addLig!(Q, l[2], 2) # put ligand at allosteric site
            for k in eachindex(slist_short)
                e_new, xyr = relaxSprings(slist_short[k], Q.pheno)
                e_old = GSE_DMS[i,j][l]
                if e_new < e_old - energyThresh
                    GSE_DMS[i,j][l] = e_new
                    structs_DMS[i,j][l] = xyr
                end
            end
            removeLig!(Q, l[1], 1) # remove ligand at active site
            removeLig!(Q, l[2], 2) # remove ligand at allosteric site
        end
        # undo the mutation
        Q = deepcopy(Qwildtype)
    end

    # test WT against new states
    for (n,l) in enumerate(ligDex)
        addLig!(Q, l[1], 1) # put ligand at active site
        addLig!(Q, l[2], 2) # put ligand at allosteric site
        for k in eachindex(slist_short)
            e_new, xyr = relaxSprings(slist_short[k], Q.pheno)
            e_old = GSE_WT[l]
            if e_new < e_old - energyThresh
                GSE_WT[l] = e_new
               # structs_DMS[i,j][l] = xyr
            end
        end
        removeLig!(Q, l[1], 1) # remove ligand at active site
        removeLig!(Q, l[2], 2) # remove ligand at allosteric site
    end

    # replace null mutations GSE's in GSE_DMS
    for (i, a) in enumerate(Q.geno)
        GSE_DMS[i,a] = GSE_WT
    end

    if saveStructs
        return GSE_WT, GSE_DMS, structs_DMS
    else
        return GSE_WT, GSE_DMS
    end
end


#function getUniqueSingleMutStructs(Q::Network,
#                                   sm_structs::Matrix{Matrix{Matrix{Float64}}})
#    
#    addLig!(Q, 1, 1) # put solvent at active site
#    addLig!(Q, 1, 2) # put solvent at allosteric site
#    structs_tmp = Matrix{Float64}[]
#    for a in sm_structs, b in a
#        if !all(b .== 0)
#            push!(structs_tmp, copy(b))
#        end
#    end
#    return removeDuplicates(Q, structs_tmp)
#end




function computeGSEAllostery(Q::Network;
                             ligDex=setLigDex("Big2"),
                             GS_settings::Tuple{Int,Int}=(10,100),
                             strains::Tuple=(-0.1,0.1),
                             energyThresh::Number=1e-5)
    # Add bonds to the back side of the network
    # and computes the change the new GSEs.

    # Requires 7x7 network
    if size(Q.A,1) != 49
        error("Cant apply computeAllostery to network of this size")
    end

    staples = [CartesianIndex(22,29),
               CartesianIndex(29,36),
               CartesianIndex(36,43),
               CartesianIndex(43,44),
               CartesianIndex(44,1),
               CartesianIndex(1,7),
               CartesianIndex(7,48),
               CartesianIndex(48,49),
               CartesianIndex(49,42),
               CartesianIndex(42,35),
               CartesianIndex(35,28)]

    numStaples = length(staples)
    Q = deepcopy(Q) # to detach from Q outside function
    removeAllLig!(Q)
    addLig!(Q, 1, 1) # solvent at active site
    addLig!(Q, 1, 2) # solvent at allosteric site
    GSE_WT = copy(Q.energies)
    Qwildtype = deepcopy(Q)
    xy_s = Q.structs[1]
    GSE_DAS = Matrix{Matrix{Float64}}(undef, numStaples, length(strains))
    slist = Matrix{Float64}[]
    for (i, bond) in enumerate(staples), (j, strain) in enumerate(strains)
        l = computeNodeNodeDist(xy_s, bond)
        Q.pheno[1][bond] = Q.pheno[1]'[bond] = 1. # set spring constant
        Q.pheno[2][bond] = Q.pheno[2]'[bond] = l * (strain + 1)  # set rest length
        Q.pheno[3][bond] = Q.pheno[3]'[bond] = 0. # set energy
        tmpGSE = zeros(size(Q.energies))
        tmpGSE[ligDex], structs = findGroundStates(Q, ligDex, GS_settings)
        GSE_DAS[i,j] = tmpGSE
        append!(slist, structs[:])
        Q = deepcopy(Qwildtype)
    end

    # only keep unique structures, compared by their energy. store them in slist_short.
    elist_short = Float64[]
    slist_short = Matrix{Float64}[]
    for i in eachindex(slist)
        e, xyr = relaxSprings(slist[i], Q.pheno)
        if !intol(e, elist_short, energyThresh) # make this a variable
            push!(elist_short, e)
            push!(slist_short, copy(xyr))
        end
    end

    # re test all of the states found from findGroundStates
    for (i, bond) in enumerate(staples), (j, strain) in enumerate(strains)
        l = computeNodeNodeDist(xy_s, bond)
        Q.pheno[1][bond] = Q.pheno[1]'[bond] = 1. 
        Q.pheno[2][bond] = Q.pheno[2]'[bond] = l * (strain * 1) 
        Q.pheno[3][bond] = Q.pheno[3]'[bond] = 0.
        for (n,l) in enumerate(ligDex)
            addLig!(Q, l[1], 1) # put ligand at active site
            addLig!(Q, l[2], 2) # put ligand at allosteric site
            for k in eachindex(slist_short)
                e_new, xyr = relaxSprings(slist_short[k], Q.pheno)
                e_old = GSE_DAS[i,j][l]
                if e_new < e_old - energyThresh
                    GSE_DAS[i,j][l] = e_new
                   # structs_DMS[i,j][l] = xyr
                end
            end
            removeLig!(Q, l[1], 1) # remove ligand at active site
            removeLig!(Q, l[2], 2) # remove ligand at allosteric site
        end
        Q = deepcopy(Qwildtype)
    end

    # test WT againts new states
    for (n,l) in enumerate(ligDex)
        addLig!(Q, l[1], 1) # put ligand at active site
        addLig!(Q, l[2], 2) # put ligand at allosteric site
        for k in eachindex(slist_short)
            e_new, xyr = relaxSprings(slist_short[k], Q.pheno)
            e_old = GSE_WT[l]
            if e_new < e_old - energyThresh
                GSE_WT[l] = e_new
               # structs_DMS[i,j][l] = xyr
            end
        end
        removeLig!(Q, l[1], 1) # remove ligand at active site
        removeLig!(Q, l[2], 2) # remove ligand at allosteric site
    end
    return GSE_WT, GSE_DAS
end



function computeGSEBondSens(Q::Network,
                           ligDex=setLigDex("Big5");
                           GS_settings::Tuple{Int, Int}=(10,100),
                           saveStructs::Bool=false,
                           energyThresh::Number=1e-5)
    
    # compute the ground state energies for a network
    # for all bond mutations
    Q = deepcopy(Q) # to detach from Q outside function
    numNodes = size(Q.A, 1)
    numBonds = Int(count(Q.A .!= 0) / 2)
    Qwildtype = deepcopy(Q)
    strains = [.2, -.2] # argument for mutateBonds
    numLigs = length(Q.ligs)
    energies = zeros(size(Q.energies))

    slist = Matrix{Float64}[]
    GSE_WT = copy(Q.energies)
    GSE_WT[ligDex], _ = findGroundStates(Q, ligDex, GS_settings)
    GSE_DBS = Matrix{Matrix{Float64}}(undef, numBonds, 2)
    bondIndices = findall(LowerTriangular(Q.A) .!=0 )

    for (i, bondIndex) in enumerate(bondIndices)
        for (j,strain) in enumerate(strains)
            mutateBond!(Q, bondIndex, strain)
            tmpGSE = zeros(size(energies))
            tmpGSE[ligDex], structs = findGroundStates(Q, ligDex, GS_settings)
            GSE_DBS[i,j] = tmpGSE
            append!(slist, structs[:]) 
            Q = deepcopy(Qwildtype)
        end
    end

    # only keep unique structures, compared by their energy. store them in slist_short.
    elist_short = Float64[]
    slist_short = Matrix{Float64}[]
    for i in eachindex(slist)
        e, xyr = relaxSprings(slist[i], Q.pheno)
        if !intol(e, elist_short, energyThresh) # make this a variable
            push!(elist_short, e)
            push!(slist_short, copy(xyr))
        end
    end


    # test all single mutants with all of the unique ground state structures.
    for (i, bondIndex) in enumerate(bondIndices)
        for (j,strain) in enumerate(strains)
             mutateBond!(Q, bondIndex, strain)
             for (n,l) in enumerate(ligDex)
                 addLig!(Q, l[1], 1) # put ligand at active site
                 addLig!(Q, l[2], 2) # put ligand at allosteric site
                 for k in eachindex(slist_short)
                     e_new, xyr = relaxSprings(slist_short[k], Q.pheno)
                     e_old = GSE_DBS[i,j][l]
                     if e_new < e_old - energyThresh
                         GSE_DBS[i,j][l] = e_new
                     end
                 end
                 removeLig!(Q, l[1], 1) # remove ligand at active site
                 removeLig!(Q, l[2], 2) # remove ligand at allosteric site
             end
             # undo the mutation
             Q = deepcopy(Qwildtype)
         end
    end

    # test WT againts new states
    for (n,l) in enumerate(ligDex)
        addLig!(Q, l[1], 1) # put ligand at active site
        addLig!(Q, l[2], 2) # put ligand at allosteric site
        for k in eachindex(slist_short)
            e_new, xyr = relaxSprings(slist_short[k], Q.pheno)
            e_old = GSE_WT[l]
            if e_new < e_old - energyThresh
                GSE_WT[l] = e_new
               # structs_DMS[i,j][l] = xyr
            end
        end
        removeLig!(Q, l[1], 1) # remove ligand at active site
        removeLig!(Q, l[2], 2) # remove ligand at allosteric site
    end
    return GSE_WT, GSE_DBS
end


function computeSingleMutsFit(GSE_DMS::Array{<:Array{Float64}}, # either Mut, Bond, Allo
                              params::Params)
    return map(x -> computeFitness(x, params), GSE_DMS )
end

function computeFitSens(GSE_WT::T,
                        GSE_DMS::Array{T}, # either Mut, Bond, Allo
                        params::Params) where T <: Array{Float64}
    return computeSingleMutsFit(GSE_DMS, params) .- computeFitness(GSE_WT, params)
end


function computeAllostery(GSE_WT::Matrix{Float64},
                          GSE_DAS::Array{<:Matrix{<:Float64}},
                          params::Params)::Float64
    return maximum(abs.(computeFitSens(GSE_WT, GSE_DAS, params)[1:end] ))
end

