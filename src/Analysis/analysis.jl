################################################
###### Single Trajectory Analysis Functions #####
#################################################


function computeStatehoodAllo(trans_mat::Matrix)::Int
    
    @assert size(trans_mat) == (4,4)
    trans_mat = Int.(trans_mat)
    t = trans_mat .& trans_mat'
    t[t^3 .!= 0] .= 1 # make the 4 node graph block diagonal
    s = sum(t)
    if s==16
        return 1
    elseif s==10 || s==8.
        return 2
    elseif s==6
        return 3
    elseif s==4.
        return 4 
    else
        error("you done goofed")
    end
end


function computeStatehood(trans_mat::Matrix)::Int
    # 
    if sum(trans_mat) == length(trans_mat)
        return 0 # single state
    else
        return 1 # multi state
    end
end


function confChangeTest(Q::Network,
                        ligDex::Vector{CartesianIndex{2}})::Int
    # method 1
    # return number of conformations network occupies
    # when binding ligands
    E_mat, confs_mat, trans_mat = swapLigs(Q, ligDex)
    statehood = computeStatehood(trans_mat)
    return statehood
end

function confChangeTest(Q::Network,
                        lig1::Ligand,
                        lig2::Ligand,
                        ligDex::Vector{CartesianIndex{2}})::Int
    # method 1
    # return number of conformations network occupies
    # when binding ligands
    Q = deepcopy(Q)
    Q.ligs[1] = lig1
    Q.ligs[2] = lig2
    E_mat, confs_mat, trans_mat = swapLigs(Q, ligDex)
    statehood = computeStatehood(trans_mat)
    return statehood
end


function confChangeTest(trans_mat::Matrix{Float64})::Int
    # using this function to find out if all of the ligand binding 
    # combinations in ligDex are in the same state.
    statehood = computeStatehood(trans_mat)
    return statehood
end


function confChangeTest2(Q::Network,
                        ligDex::Vector{CartesianIndex{2}})::Int
    # Use this function for finding out if the network
    # is 1,2,3 or 4 state allostery.
    E_mat, confs_mat, trans_mat = swapLigs(Q, ligDex)
    statehood = computeStatehoodAllo(trans_mat)
    return statehood
end


function computeLeastStepMatrix(adjacencyMatrix_in)
    # compute a matrix where each element is the least number of steps between 
    # the i,jth nodes.
    
    adjacencyMatrix = convert(Array{Float64,2},adjacencyMatrix_in .!= 0)

    numNodes = size(adjacencyMatrix,1)
    adjacencyMatrix = adjacencyMatrix + eye(numNodes)
    leastStepMatrix = copy(adjacencyMatrix)
    n=1
    while any(leastStepMatrix .== 0)
        n += 1
        Alow = adjacencyMatrix^(n-1)
        low = Alow .!= 0
        Ahigh = Alow * adjacencyMatrix
        high = Ahigh .!= 0
        typeof(high)
        indicies = high .& .!low
        leastStepMatrix[indicies] .= n
    end
    leastStepMatrix[eye(numNodes).>0] .= 0
    return leastStepMatrix
end



function computeStepsFromAct(leastStepMatrix; actNodes=[2,3,4])
    # computes the least number of steps from each node to 
    # a the defined active site nodes {2,3,4}.
    stepsVec = vec(minimum(leastStepMatrix[:,actNodes], dims=2))
    return stepsVec
end



function computeBondStrains(xy1, xy2, adjacencyMatrix)
    # xy1 is the solvent bound network positions
	# xy2 is the ligand bound network positions
    trueAdjMat = adjacencyMatrix .!= 0
    DM1 = computeDistanceMatrix(xy1, xy1)
	DM2 = computeDistanceMatrix(xy2, xy2)
    strainMatrix = (DM2 .- DM1) ./ ((DM1 .+ DM2) ./ 2) 
    strainMatrix .= trueAdjMat .* strainMatrix
    return strainMatrix
end


function locateConfChange(Q::Network,
                          ligDex::Vector{CartesianIndex{2}};
                          numBondThresh::Int=4 # top most strain bonds
                         )::Float64

    if length(ligDex) != 2
        error("Error: ligDex does not have the right size")
    end
    steps = zeros(numBondThresh)
    SM, RLM, EM = deepcopy(Q.pheno)
    E_mat, confs_mat, trans_mat = swapLigs(Q, ligDex)
   
   # test is there is a conformational change
    statehood = confChangeTest(trans_mat)
    if statehood == 0 
        return -1.
    end
    
    # find top strained bonds
    tmp1 = LowerTriangular(abs.(computeBondStrains(confs_mat[1], confs_mat[2], Q.A)))[:]
    bonds = CartesianIndices(size(Q.A))[sortperm(tmp1, rev=true)[1:numBondThresh]]
    # compute steps from active site
    LSM = computeLeastStepMatrix(Q.A)
    sfas = computeStepsFromAct(LSM; actNodes=collect(Q.sites[1][1][1]:Q.sites[1][1][2]))
    for i in 1:numBondThresh
        steps[i] = minimum(sfas[collect(Tuple(bonds[i]))])
    end
    return maximum(steps)
end


function extractEnergyTrace(Qs::Vector{Network})::Matrix{Float64}
    # extracts data from a single evolutionary trajectory
    
    numMuts = length(Qs)
    energyArray = zeros(size(Qs[1].energies)..., numMuts)
    for i in 1:numMuts
        energyArray[:,:,i] = Qs[i].energies
    end
    return energyArray
end


function extractEnergyTrace(energiesList::Vector{Matrix{Float64}})
    # extracts data from a single evolutionary trajectory
    
    numMuts = length(energiesList)
    energyArray = zeros(size(energiesList[1])..., numMuts)
    for i in 1:numMuts
        energyArray[:,:,i] = energiesList[i]
    end
    return energyArray
end

function computeNodeNodeDist(xy::Array{Float64,2},
                             i::Int64,
                             j::Int64)::Float64
    dx = xy[i,1] - xy[j,1]
    dy = xy[i,2] - xy[j,2]
    return sqrt(dx^2 + dy^2)
end

function computeNodeNodeDist(xy::Array{Float64,2}, ij::CartesianIndex)::Float64
    return computeNodeNodeDist(xy, ij[1], ij[2])
end


function computeEnergyLigSlice(Q::Network,
                               ligLengths::AbstractVector{<:Number};
                               numTrials::Int=1000)

    N = length(ligLengths)
    energies1 = zeros(N)
    energies2 = zeros(N)

   # structs = Matrix{Matrix{Float64}}(undef, N,N)
    xy = Q.structs[1] # 
    ac=Q.sites[1][1] # this only uses the first bond of the sites
    al=Q.sites[2][1]
    SM, RLM, EM = Q.pheno
    SM[ac] = SM'[ac] = 1.0
    SM[al] = SM'[al] = 1.0
    EM[ac] = EM'[ac] = 0.0
    EM[al] = EM'[al] = 0.0

    # First with the solvent at the allosteric site
    RLM[al] = RLM'[al] = Q.ligs[1].bonds[1].l
    for (i,l) in enumerate(ligLengths)
        RLM[ac] = RLM'[ac] = l # manually define ligand length
        energies1[i], _ = findGroundState(xy, SM, RLM, EM, numTrials=numTrials)
    end

    # Second with the ligand at the allosteric site
    RLM[al] = RLM'[al] = Q.ligs[2].bonds[1].l
    for (i,l) in enumerate(ligLengths)
        RLM[ac] = RLM'[ac] = l # manually define ligand length
        energies2[i], _ = findGroundState(xy, SM, RLM, EM, numTrials=numTrials)
    end

    return energies1, energies2
end


function computeLigScape(Q::Network,
                        ligLengths::AbstractVector{<:Number};
                        GS_settings::Tuple{Int, Int} = (10, 100),
                        energyThresh::Number=1e-5)

    N = length(ligLengths)
    energies = zeros(N, N)
    structs = Matrix{Matrix{Float64}}(undef, N,N)
    xy = buildxy0()
    ac=Q.sites[1][1] # this only uses the first bond of the sites
    al=Q.sites[2][1]
    SM, RLM, EM = Q.pheno
    SM[ac] = SM'[ac] = 1.0
    SM[al] = SM'[al] = 1.0
    for (i1,l1) in enumerate(ligLengths), (i2,l2) in enumerate(ligLengths)
        RLM[ac] = RLM'[ac] = l1 # manually define ligand length
        RLM[al] = RLM'[al] = l2 # manually define ligand length
        energies[i1,i2], structs[i1,i2] = findGroundState(xy, SM, RLM, EM, GS_settings)
    end
    
    # only keep unique structures, compared by their energy. store them in slist_short.
    elist_short = Float64[]
    slist_short = Matrix{Float64}[]
    xy_list = [structs[:]; Q.structs[[1,2,4,5]]]
    removeAllLig!(Q)
    for i in eachindex(xy_list)
        e, xyr = relaxSprings(xy_list[i], Q.pheno)
        if !intol(e, elist_short, energyThresh) # make this a variable
            push!(elist_short, e)
            push!(slist_short, copy(xyr))
        end
    end

    # test all of ligand space against all of the unique structs.
    SM[ac] = SM'[ac] = 1.0
    SM[al] = SM'[al] = 1.0
    for (i1,l1) in enumerate(ligLengths), (i2,l2) in enumerate(ligLengths)
        RLM[ac] = RLM'[ac] = l1 # manually define ligand length
        RLM[al] = RLM'[al] = l2 # manually define ligand length
        E = energies[i1,i2]
        for xy in slist_short
            e, xy_r = relaxSprings(xy, Q.pheno)
            if e < E - energyThresh
                E = e
                energies[i1,i2] = e
                structs[i1,i2] = xy_r
            end
        end
    end
    return energies, structs
end


function labelLigScape(Q::Network, structs::Array{Matrix{Float64}})
    # find the regions in ligand space that are in the same state.

    function removeLigsAndRelax(Q,structs)
        Q = deepcopy(Q)
        structs = deepcopy(structs)
        energies = zeros(size(structs))
        removeAllLig!(Q)
        for i in eachindex(structs)
            energies[i], structs[i] = relaxSprings(structs[i], Q.pheno)
        end
        return energies
    end
    
    function labelLigScape(v::Array)
        uniqs = [v[1]]
        for i in 2:length(v)
            e = v[i]
            if !intol(e, uniqs, 1e-5)
                push!(uniqs, e)
            end
        end
        labes = zeros(Int, size(v))
        for i in 1:length(v)
            _, labes[i] = findmin(abs.(uniqs .- v[i]))
        end
        return labes
    end

    energies = removeLigsAndRelax(Q,structs)
    labes = labelLigScape(energies)
    return labes
end


function outlineLigScape(stateLabels::Matrix{Int}, ligLengths::AbstractVector)
    # Draw boundaries of regions in ligand space of same state.  
    stateLabels = reverse(stateLabels, dims=1)
    verticals = diff(stateLabels, dims=2) # vertical edges
    horizontals = diff(stateLabels, dims=1) # horizonal edges
    xs = []
    ys = []
    n = length(ligLengths)
    l0 = ligLengths[1]
    l1 = ligLengths[end]
    dl = abs(l1 - l0) / (n-1)
    
    # collect verticals segments
    for i in 1:size(verticals, 1), j in 1:size(verticals,2)
        if verticals[i,j] != 0
            x = [j*dl, j*dl] .+ l0 .- dl/2
            y = l1 .- [(i-1)*dl, i*dl] .+ dl/2
            push!(xs,x)
            push!(ys,y)
        end
    end
    
    # collect verticals segments
    for i in 1:size(horizontals, 1), j in 1:size(horizontals,2)
        if horizontals[i,j] != 0
            x = [(j-1)*dl, j*dl] .+ l0 .- dl/2
            y = l1 .- [i*dl, i*dl] .+ dl/2
            push!(xs,x)
            push!(ys,y)
        end
    end
    return xs, ys
end

gaussian(x,μ,σ) = (1/(σ*sqrt(2*pi))) * exp(-0.5*((x-μ)/σ)^2)  

#function doubleGaussian(x, μ1, μ2, σ1, σ2, w1, w2)
#    return w1*gaussian(x, μ1, σ1) + w2*gaussian(x, μ2, σ2)
#end
#
#function computeBimodality(dgparams::Tuple)
#    μ1, μ2, σ1, σ2, w1, w2 = dgparams
#    x1, x2 = μ1, μ2
#    y1 = doubleGaussian(x1, dgparams...)
#    y2 = doubleGaussian(x2, dgparams...)
#    y_mid = (y1+y2)/2
#    y_dg = doubleGaussian((x1+x2)/2, dgparams...)
#    return y_mid - y_dg
#end
    
    
function getDeformation(xy1::Matrix{Float64},
                        xy2::Matrix{Float64},
                        Q::Network)
    
    removeAllLig!(Q)
    _, xy1_relaxed = relaxSprings(xy1, Q.pheno, gradTol = 1e-10)
   # _, xy2 = relaxSprings(xy2, pheno, gradTol = 1e-10)

    H = computeHessian(xy1_relaxed, Q.pheno)
    E = eigen(H)
    x1 = xy1'[:]
    x2 = xy2'[:]
    dx = x1 .- x2
    @assert sum(abs.(E.values) .< 1e-4) == 3
    dx .-= (dx ⋅ E.vectors[:,1]) .* E.vectors[:,1]  
    dx .-= (dx ⋅ E.vectors[:,2]) .* E.vectors[:,2] 
    dx .-= (dx ⋅ E.vectors[:,3]) .* E.vectors[:,3] 
    return Array(reshape(dx, 2, :)')
end




function computeOverlapsStiffnesses(Q::Network)
    orientNetwork!(Q)
    removeAllLig!(Q)
    xy1 = Q.structs[1]
    xy2 = Q.structs[5]
    _, xy1_relaxed = relaxSprings(xy1, Q.pheno, gradTol = 1e-10)
    H = computeHessian(xy2x(xy1_relaxed), Q.pheno[1], Q.pheno[2])
    E = eigen(H)
    @assert sum(abs.(E.values) .< 1e-4) == 3
    x1 = xy1'[:]
    x2 = xy2'[:]
    dx = x1 .- x2
    dx .-= (dx ⋅ E.vectors[:,1]) .* E.vectors[:,1]  
    dx .-= (dx ⋅ E.vectors[:,2]) .* E.vectors[:,2] 
    dx .-= (dx ⋅ E.vectors[:,3]) .* E.vectors[:,3]
    normalize!(dx)
    return E.vectors' * dx, E.values
end

function computeOverlaps(Q::Network)
    overlaps, stiffnesses = computeOverlapsStiffnesses(Q)
    return overlaps
end



function computeCCSize(dxy::Matrix{Float64})
    # returns the magnitude of the conformation change.
    return norm(dxy'[:])
end


function measureInterFrac(Q::Network,k::Int)
    # this function applies to node competition experiments
    # this function count the fraction of bonds that have a 
    # rest length defined by the upper left block of the 
    # rest length table.
    # k is the size of the disordered block.
    c = 0
    numBonds = Int(sum(Q.A) / 2)
    for i in eachindex(Q.geno)
        type1 = Q.geno[i]
        nodes = findall(Q.A[:,i] .!= 0)
        for j in eachindex(nodes)
            type2 = Q.geno[nodes[j]]
            if !(type1 >= k+1 || type2 >= k+1)
                c += 1
            end
        end
    end
    c = Int(c / 2) 
    measuredInterFrac = c / numBonds
    return measuredInterFrac
end


function countSoftBonds(Q::Network;
                        k_soft = 0.01)
    return Int(count(Q.pheno[1] .== k_soft) / 2)
end



function computeFracOfSites(df::Matrix{Float64}, sign::Function, thresh::Number)
    numMuts = size(df, 2)
    df = mean(df, dims=2) * (numMuts / (numMuts - 1))
    return count(sign.(df,thresh)) / size(df,1)
end


function measureBarrier(Q::Network; withLigands::Bool=true)

    SM, RLM, EM = Q.pheno
    xy1, xy2 = Q.structs[[1,5]]
    dxy = xy2 .- xy1
   
    if withLigands

        # add solvents
        addLig!(Q, 1, 1); addLig!(Q, 1, 2)
        E_xy1 = computeEnergy(xy1, SM, RLM, EM)
        E_xy2 = computeEnergy(xy2, SM, RLM, EM)
        E_mid = computeEnergy(xy1 .+ 0.5 .* dxy , SM, RLM, EM)
        barrierHeight1 = E_mid - mean([E_xy1,E_xy2])

        # add ligands
        addLig!(Q, 2, 1); addLig!(Q, 2, 2)
        E_xy1 = computeEnergy(xy1, SM, RLM, EM)
        E_xy2 = computeEnergy(xy2, SM, RLM, EM)
        E_mid = computeEnergy(xy1 .+ 0.5 .* dxy , SM, RLM, EM)
        barrierHeight2 = E_mid - mean([E_xy1,E_xy2])
        bh =  ( barrierHeight1 + barrierHeight2 ) / 2

    elseif ~withLigands
        # add solvents
        removeAllLig!(Q)
        n = 501
        x = LinRange(-0.1, 1.1, n)
        energies = zeros(n)
        for i in 1:n
            energies[i] = computeEnergy(xy1 .+ x[i] .* dxy, SM, RLM, EM)
        end

        roots = findRoots(diff(energies)) 
        curve = diff(diff(energies))[roots]
        extrema = roots .+ 1

        if length(roots) <= 1
            bh = zero(Float64)
        elseif length(roots) >= 2
            minima = extrema[findall( curve .> 0)]
            maxima = extrema[findall( curve .< 0)]
            bh  = maximum(energies[maxima]) - maximum(energies[minima])
        end
    end
    return bh
end


function findRoots(v::Vector)
    #find the index i where v_i~0
    roots = Int[]
    for i in 1:length(v)-1
        v_left = v[i]
        v_right = v[i+1]
        if v_left*v_right < 0
            push!(roots, i)
        end
     end
     return roots
end

function computeEnergyAlongCoor(Q::Network)

    SM, RLM, EM = Q.pheno
    xy1 = Q.structs[1]
    xy2 = Q.structs[5]
    dxy = xy2 .- xy1
    N = 100
    x = LinRange(1,1.3, N)

    removeAllLig!(Q)
    energies = zeros(N)
    for (i,a) in enumerate(x)
        xy = xy1 .+ a .* dxy 
        energies[i] = computeEnergy(xy, SM, RLM, EM)
    end
    return energies
end


function computeRuggedness(Q::Network;
                           N_reps::Int=1000,
                           pertMag::Number=4,
                           tol=1e-5)

    removeAllLig!(Q) 
    SM, RLM, EM = Q.pheno
    x0 = xy2x(buildxy0())
    energies = Float64[]
    for i in 1:N_reps
        x_perturbed = applyRandStrain(x0, pertMag)
        E_relaxed, xy_relaxed = relaxSprings(x_perturbed, SM, RLM, EM, gradTol=1e-10)
        if !intol(E_relaxed, energies, tol) 
            push!(energies, E_relaxed)
        end
    end
    return energies

end


