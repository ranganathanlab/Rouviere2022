
function buildNetwork(S::Dict)

    spacing = 1
    wid = S["wid"]
    len = S["len"]
    rando = S["rando"]

    # Build disordered hexagoanl lattice
    # return positions and adjacency matrix

    # for proper network orientation, Temp fix
    len == 7 ? a = -1 : error("len can only be 7 right now")

    numNodes = wid*len
    xvals = zeros(numNodes,1)
    yvals = copy(xvals)
    pxvals = copy(xvals)
    pyvals = copy(xvals)

    n = 0
    # construct lattice.    
    for i = 1:wid, j = 1:len
        n += 1
        xvals[n] = a*rem(j,2)*spacing/2 + (i*spacing) + 2*rando*(rand()-0.5)
        yvals[n] = j*spacing*(sin(acos(0.5))) + 2*rando*(rand()-0.5)
        pxvals[n] = a*rem(j,2)*spacing/2 + (i*spacing)
        pyvals[n] = j*spacing*(sin(acos(0.5)))
    end

    # the special nodes
    i=8; j=3
    xvals[1] = a*rem(j,2)*spacing/2 + (i*spacing) + 2*rando*(rand()-0.5)
    yvals[1] = j*spacing*(sin(acos(0.5))) + 2*rando*(rand()-0.5)
    pxvals[1] = a*rem(j,2)*spacing/2 + (i*spacing)
    pyvals[1] = j*spacing*(sin(acos(0.5)))

    # the special nodes
    i=8; j=5
    xvals[7] = a*rem(j,2)*spacing/2 + (i*spacing) + 2*rando*(rand()-0.5)
    yvals[7] = j*spacing*(sin(acos(0.5))) + 2*rando*(rand()-0.5)
    pxvals[7] = a*rem(j,2)*spacing/2 + (i*spacing)
    pyvals[7] = j*spacing*(sin(acos(0.5)))

    # build Adjacency Matrix
    A = zeros(numNodes,numNodes);
    for i in 2:numNodes
        for j in 1:i-1
            dx = pxvals[i] - pxvals[j]
            dy = pyvals[i] - pyvals[j]
            dist = sqrt(dx^2 + dy^2)
            if dist < 1.1
                A[i,j] = A[j,i] = 1.
            end
        end
    end

    #xy = zeros(length(xvals),2)
    #xy[:,1] = xvals
    #xy[:,2] = yvals
    xy_hexLattice = [pxvals pyvals]
    xy = [xvals yvals]
    NRLM = (computeDistanceMatrix(xy) .- computeDistanceMatrix(xy_hexLattice)) .* A
    return xy, A, NRLM
end


function buildxy0()
    # This is a temporary function
    spacing = 1
    wid = 7
    len = 7
    numNodes = wid*len
    xvals = zeros(numNodes,1)
    yvals = copy(xvals)
    a = -1
    n = 0
    # construct lattice.    
    for i = 1:wid, j = 1:len
        n += 1
        xvals[n] = a*rem(j,2)*spacing/2 + (i*spacing)
        yvals[n] = j*spacing*(sin(acos(0.5)))
    end
    # the special nodes
    i=8; j=3
    xvals[1] = a*rem(j,2)*spacing/2 + (i*spacing)
    yvals[1] = j*spacing*(sin(acos(0.5)))
    # the special nodes
    i=8; j=5
    xvals[7] = a*rem(j,2)*spacing/2 + (i*spacing)
    yvals[7] = j*spacing*(sin(acos(0.5)))
    return [xvals yvals]
end

function assignGenotype(numTypes::Int, numNodes::Int)
    # initialize a seqence
    genotype = rand(1:numTypes, numNodes)
    return genotype
end


function buildGeno2PhenoTable(numTypes::Int,
                              num_soft::Int,
                              k_soft::Number,
                              restLengthRange::Float64,
                              E_soft::Float64;
                              isRestLengthBlocky::Bool=false,
                              sizeRestLengthBlock::Int=5
                             )::Array{Float64,3}
    # build a table that stores the genotype to phenotype 
    # map. spring constant, rest length, atractive energy
    # data for each bond.

    rng = MersenneTwister(999999999); # dont use this seed after
    
    typeTable = zeros(numTypes, numTypes, 3)

    # Spring constants: first slice
    tmp = ones(numTypes, numTypes)
    inds = findall(LowerTriangular(tmp) .> 0)
    tmp[inds[ randperm(length(inds))[1:num_soft] ]] .= k_soft
    #tmp[rand(rng, numTypes, numTypes) .< frac_soft] .= k_soft
    typeTable[:,:,1] .= symmetric(tmp)

    # Rest Length factor: second slice
    typeTable[:,:,2] .= symmetric( (restLengthRange*
                        #rand(rng, numTypes, numTypes)) .+ 1 .- restLengthRange/2)
                        rand(numTypes, numTypes)) .+ 1 .- restLengthRange/2)

    # Energy of interaction: third slice
    typeTable[:,:,3] = E_soft * Float64.(typeTable[:,:,1] .!= 1)
    #typeTable[:,:,3] = std_E .* symmetric(randn(rng, numTypes, numTypes)) .+ mean_E


    # for genotype competition table
    if isRestLengthBlocky
        # fix the rest length
        tmp = ones(numTypes, numTypes)
        tmp[1:sizeRestLengthBlock,1:sizeRestLengthBlock] = typeTable[1:sizeRestLengthBlock,1:sizeRestLengthBlock,2]
        typeTable[:,:,2] = tmp

        # fix the spring constants
        tmp = ones(numTypes, numTypes)
        tmp2 = ones(numTypes-sizeRestLengthBlock, numTypes - sizeRestLengthBlock)
        tmp3 = ones(numTypes-sizeRestLengthBlock, numTypes - sizeRestLengthBlock)
        tmp4 = ones(numTypes-sizeRestLengthBlock, numTypes - sizeRestLengthBlock)
        inds = findall(LowerTriangular(tmp2) .> 0)
        
        tmp2[inds[ randperm(length(inds))[1:num_soft] ]] .= k_soft
        tmp3[inds[ randperm(length(inds))[1:num_soft] ]] .= k_soft
        tmp4[randperm(length(tmp4))[1:num_soft] ] .= k_soft
        tmp[sizeRestLengthBlock+1:end, sizeRestLengthBlock+1:end] = symmetric(tmp2)
        tmp[1:sizeRestLengthBlock, 1:sizeRestLengthBlock] = symmetric(tmp3)
        tmp[1:sizeRestLengthBlock, sizeRestLengthBlock+1:end] = tmp4
        tmp[sizeRestLengthBlock+1:end, 1:sizeRestLengthBlock] = tmp4'
        typeTable[:,:,1] .= tmp
        
    end
    return typeTable
end



function geno2Pheno(typeTable, genotype, adjacencyMatrix, NRLM)
    # Maps the genotype to a phenotype that is stored in phenoMatrix

    numNodes = length(genotype)
    numTypes = size(typeTable,1)
    z = zeros(numNodes,numNodes)
    phenotype = [copy(z), copy(z), copy(z)]
    for i in 1:numNodes
        for j in 1:numNodes
            typei = genotype[i]
            typej = genotype[j]
            # spring matrix
            phenotype[1][i,j] = adjacencyMatrix[i,j] * typeTable[typei,typej,1]
            # rest length matrix
            #phenotype[2][i,j] = adjacencyMatrix[i,j] * typeTable[typei,typej,2]
            phenotype[2][i,j] = adjacencyMatrix[i,j] * typeTable[typei,typej,2] + NRLM[i,j]
            # energy of interaction matrix
            phenotype[3][i,j] = adjacencyMatrix[i,j] * typeTable[typei,typej,3]
        end
    end
    return phenotype
end


# the new ones
function buildNetworkFromSettings(S::Dict)::Network
    # Build Network Object
    numNodes = S["len"] * S["wid"]
    xy0, A, NRLM  = buildNetwork(S)
    typeTable = buildGeno2PhenoTable(S["numTypes"], S["num_soft"], S["k_soft"], S["restLengthRange"], S["E_soft"],
                                     isRestLengthBlocky=S["isRestLengthBlocky"],
                                     sizeRestLengthBlock=S["sizeRestLengthBlock"])

    genotype = assignGenotype(S["numTypes"], numNodes)
    phenotype =  geno2Pheno(typeTable, genotype, A, NRLM)
    numLigs = length(S["ligs"])
    energies = zeros(numLigs, numLigs)
    structs = Matrix{Matrix{Float64}}(undef, numLigs, numLigs)
    [ structs[i] = deepcopy(xy0) for i in eachindex(structs)  ]
    Q = Network(A, NRLM, genotype, phenotype, typeTable, energies, structs,
                S["ligs"], S["sites"], S["unfoldedEnergy"], S["suppFitnessParams"])
    return Q
end


function buildNetworkFromSeq(A::Matrix{Float64},
                             geno::Vector{Int64}, 
                             typeTable::Array{Float64,3}, 
                             ligs::Vector{Ligand},
                             sites::Vector{Vector{CartesianIndex{2}}},
                             unfoldedEnergy::Float64,
                             energies=nothing, 
                             structs=nothing)
    # with this function I only need to save the sequence
    # and I can get back the rest.
    numSites = length(sites)
    numLigs = length(ligs)
    numNodes = length(geno)
    numDims = 2
    if energies == nothing
        energies = zeros(numLigs,numLigs)
    end
    if structs == nothing
        structs =[zeros(numNodes, numDims) for i in 1:numLigs, j in 1:numLigs]
    end
    pheno = geno2pheno(typeTable, geno, A)
    return Network(A, geno, pheno, typeTable, energies, structs, ligs, sites, unfoldedEnergy)  
end


function buildNetworkFromSeq(geno::Vector{Int64},
                             Q::Network,
                             energies=nothing, 
                             structs=nothing)
    numSites = length(Q.sites)
    numLigs = length(Q.ligs)
    numNodes = length(geno)
    numDims = 2

    if energies == nothing
        energies = zeros(numLigs,numLigs)
    end
    if structs == nothing
        structs =[zeros(numNodes, numDims) for i in 1:numLigs, j in 1:numLigs]
    end

    pheno = geno2pheno(Q.typeTable, geno, Q.A)
    return buildNetworkFromSeq(A, geno, Q.typeTable, Q.ligs,
            Q.sites, Q.unfoldedEnergy, energies=energies, structs=structs)  
end


function buildNetworkFromSeq(L::LilNetwork,
                             Q::Network)
    buildNetworkFromSeq(L.geno, Q, energies=L.energies, structs=L.structs)
end

function makeLilNetwork(Q::Network)
    return LilNetwork(Q.geno, Q.energies, Q.structs)
end



function convertNetwork2LiteNetwork(Q::Network, ligDex::Vector{CartesianIndex{2}})::LiteNetwork
    # convert a Network to a LiteNetwork to save space when saving.
    A = BitMatrix(Q.A)
    NRLM = convertSymMat2Vec(Q.NRLM)
    #pheno = [convertSymMat2Vec(Q.pheno[i]) for i in eachindex(Q.pheno)]
    energies = Q.energies[ligDex]
    structs = Q.structs[ligDex]
    preligs = ligs2preligs(Q.ligs)
    return LiteNetwork(A, NRLM, Q.geno, Q.typeTable, energies,
        structs, preligs, Q.sites, Q.unfoldedEnergy, Q.suppFitnessParams, ligDex)
end

function convertLiteNetwork2Network(LQ::LiteNetwork)::Network
    # convert a LiteNetwork to a Network.
    
    numLigs = length(LQ.preligs)
    
    A = Matrix{Float64}(LQ.A)
    NRLM = convertVec2SymMat(A, LQ.NRLM)
    pheno = geno2Pheno(LQ.typeTable, LQ.geno, A, NRLM)
    energies = zeros(numLigs, numLigs) 
    energies[LQ.ligDex] = LQ.energies
    structs = Matrix{Matrix{Float64}}(undef, numLigs, numLigs)
    [ structs[i] = deepcopy(buildxy0()) for i in eachindex(structs)  ]
    structs[LQ.ligDex] .= LQ.structs
    ligs = preligs2ligs(LQ.preligs)
    return Network(A, NRLM, LQ.geno, pheno, LQ.typeTable, energies, structs,
            ligs, LQ.sites, LQ.unfoldedEnergy, LQ.suppFitnessParams)
    
end






                             
