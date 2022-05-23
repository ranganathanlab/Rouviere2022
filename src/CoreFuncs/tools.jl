##########################################################
####  General tools for model ############################
##########################################################


function eye(n)    
    # like the eye() in matlab, uses LinearAlgebra
    #Matrix{Float64}(I, n, n)
    diagm(0=>ones(n))
end

function acceptMC(delta_E::Number, T::Number) #
    # Metropolis monte carlo

    if delta_E <= 0
        return true
    elseif delta_E > 0
        p = exp( -delta_E / T )
        u = rand()
        if u <= p
            return true
        else
            return false
        end
    end
end

function symmetric!(mat)
    # takes a matrixs and copies lower triangle
    # to the upper triangle therebye making it
    # symmetric. Diagonal remains same.

    for i in 2:size(mat,1)
        for j in 1:i-1
            mat[j,i] = mat[i,j]
        end
    end
    return nothing
end



function symmetric(mat)
    # takes a matrixs and copies lower triangle
    # to the upper triangle therebye making it
    # symmetric. Diagonal remains same.

    MAT = copy(mat)
    symmetric!(MAT)
    return MAT
end

function linspace(srt, stp, n)
    # like matlab's linspace
    return collect(range(srt,stop=stp,length=n))
end


function movFromPics()
    cmd = `convert -delay 6 -quality 80 '*.png' ../movie.mp4`
    run(cmd)
    return nothing
end


function computeProgress(current, total)
    per = Float64(floor(100*current/total))
    message = @sprintf("%4.0f", per)
    return message
end


function displayProgress(current, total)
    message = computeProgress(current, total)
    print(" $message% done\r")
end


function makeIfNot(dir)
    if !isdir(dir)
        mkdir(dir)
    end; return nothing
end


function l2a(list)
    # convert a list of lists to a dim2 Array
    
    N = length(list)
    # check that each element list is of same length
    M = length(list[1])
    for i in 2:N
        test_M = length(list[i])
        if test_M != M
            error("Not all element list are of the same length")
        end
    end

    A = Array{typeof(list[1][1]),2}(undef, N, M)
    for i in 1:N, j in 1:M
        A[i,j] = list[i][j]
    end
    return A
end


function issquare(M::AbstractMatrix)::Bool
    if size(M,1) == size(M,2)
        return true
    else
        return false
    end
end


function k2ij(k::Int, N::Int)
    # return the i,j indicies of a lower-triangular matrix at poistion k.
    # N is the sidelength of the square matrix.
    @assert k <= (N^2 - N) / 2
    notFinished = true
    j = 0
    summ = 0
    i = 0
    while notFinished
        j += 1
        summ += N - j
        if k <= summ
            i = N - summ + k
            notFinished = false
        end
    end
    return i,j
end



function vec2LowerTriMat!(M::AbstractMatrix{T},
                          v::Vector{T}
                          ) where T<:Number

    # copy contents of vector to the lower triangle
    # of M not including diagonal.
    @assert issquare(M)
    
    # make sure that v is of the correct length
    n = size(M,1)
    @assert k2ij(length(v),n) == (n, n-1)
    
    for k in 1:length(v)
        i,j = k2ij(k, n)
        M[i,j] = v[k]
    end
    return nothing
end


function vec2UpperLowerTriMat!(M::AbstractMatrix{T},
                               v::Vector{T}
                               ) where T<:Number

    # copy contents of vector to the upper 
    # and lower triangle of M not including diagonal.

    vec2LowerTriMat!(M,v)
    vec2LowerTriMat!(M',v)
    return nothing
end


function randp(P::Vector{AbstractFloat})::Int
    # return an index of the argument vector
    # with probability defined by the elements
    # of the vector.

    @assert sum(P) == 1.0
    cumm=0.0
    r = rand()
    for (i,p) in enumerate(P)
        cumm += p
        if r <= cumm
            return i
        end
    end
end



function intol(x::Number, v::Vector{<:Number}, tol::Number)
    # returns true if there is a number in v that is 
    # within +/- tol of x.
    for y in v
        y_plus = y + tol
        y_minus = y - tol
        if y_minus <= x <= y_plus
            return true
        end
    end
    return false
end





function xy2x(xy::Matrix{Float64})
    return xy'[:]
end

function x2xy(x::Vector{Float64})
    return collect(reshape(x, 2, :)')
end

function computeDistanceMatrix(xy1::Matrix{Float64},
                               xy2::Matrix{Float64}
                              )::Matrix{Float64}
    # this function compute the distance matrixs between points of xy1 and
    # points of xy2. the (i,j) entry of the matrix is the distance from the
    # ith point of xy1 to the jth point of xy2.

    N = size(xy1,1)
    distMat = zeros(N,N)
    @inbounds for i in 1:N
        xi = xy1[i,1]
        yi = xy1[i,2]
        @simd for j in 1:N
            dx = xi - xy2[j,1]
            dy = yi - xy2[j,2]
            distMat[i,j] = sqrt(dx^2 + dy^2)
        end
    end
    return distMat
end


function computeDistanceMatrix(xy::Matrix{Float64})::Matrix{Float64}
    return computeDistanceMatrix(xy,xy)
end


function computeNodeNodeDistance(xy::Matrix{Float64},
                                 node1::Int,
                                 node2::Int)
    x1, y1 = xy[node1, :]
    x2, y2 = xy[node2, :]
    return sqrt( (x2-x1)^2  + (y2-y1)^2)
end

function computeConfDiff(xy1::Matrix{Float64},
                         xy2::Matrix{Float64})::Float64
    # this function compute the distance matrixs between points of xy1 and
    # points of xy2. the (i,j) entry of the matrix is the distance from the
    # ith point of xy1 to the jth point of xy2.

    N = size(xy1,1)
    diff = 0.0 
    @inbounds for i in 1:N
        xi1 = xy1[i,1]
        yi1 = xy1[i,2]
        xi2 = xy2[i,1]
        yi2 = xy2[i,2]
        @simd for j in 1:N
            dx1 = xi1 - xy1[j,1]
            dy1 = yi1 - xy1[j,2]
            dx2 = xi2 - xy2[j,1]
            dy2 = yi2 - xy2[j,2]
            d = abs( sqrt(dx1^2 + dy1^2) - sqrt(dx2^2 + dy2^2))
            diff = ifelse(d>diff, d, diff)
        end
    end
    return diff
end


function sameConf(xy1::Matrix{Float64},
                  xy2::Matrix{Float64})::Bool
    # answers the question are these two structures
    # the same conformation
    if  computeConfDiff(xy1, xy2) < .1 # theshold chosen from analysis
        return true
    else
        return false
    end
end


function findLowEIndices(E_table::Matrix{Float64})
    # find the vertical indices of the lowest 
    # val of each column of E_table.
    indices = [e[1] for e in argmin(E_table, dims=1)]
    return indices
end

############################################


function center!(x::Vector{Float64})
    # put the centroid of the array
    n = Int(length(x)/2)

    avg_x = 0.0
    avg_y = 0.0
    for i in 1:n
        avg_x += x[2n-1]
        avg_y += x[2n]
    end

    avg_x /= n
    avg_y /= n

    for i in 1:n
        x[2i-1] -= avg_x
        x[2i] -= avg_y
    end
    return nothing
end


function center!(xy::Matrix{Float64})
    # put the centroid of the array
    xx = view(xy,:,1)
    yy = view(xy,:,2)
    xx .-= mean(xx)
    yy .-= mean(yy)
    return nothing
end


function center(x::Vector{Float64})::Vector{Float64}
    x = copy(x)
    center!(x)
    return x
end


function rotate!(xy::Matrix{Float64}, θ::Float64)
    # rotate network counter clockwise θ radians
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    for i in 1:size(xy,1)
      #  xy[i,:] = xy[i,:]' * R
        xy[i,:] =  R * xy[i,:]
    end
    return nothing
end


function measureAngle(xy::Matrix{Float64}, ASB::CartesianIndex{2})::Float64
    # measure how rotated the network is,
    # Must be done after network has been centered
    xy_avg = vec(mean([xy[ASB[1],:] xy[ASB[2],:]]', dims=1))
    r = norm(xy_avg)
    θ = pi - acos(xy_avg[1]/r)
    if xy_avg[2] >= 0
        return θ
    else
        return -θ
    end
end


function orient!(xy::Matrix{Float64}, ASB::CartesianIndex{2})
    center!(xy)
    θ = measureAngle(xy, ASB)
    rotate!(xy, θ)
    return nothing
end

function orientNetwork!(Q::Network)
    for i in 1:length(Q.structs)
        orient!(Q.structs[i], Q.sites[1][1])
    end
    return nothing
end



function removeTranslation!(dx::Vector{Float64})
    # removes translations from the perturbation dxy
    center!(dx)
    return nothing
end

function removeRotation!(dx::Vector{Float64}, x::Vector{Float64})
    # removes rotaions from the perturbation dxy
    
    n = Int(length(x)/2)
    xx = x[1:2:end]
    yy = x[2:2:end]
    avgX = mean(xx)
    avgY = mean(yy)

    rotationVec = similar(x)
    for i in 1:n
        rotationVec[2i-1] = x[2i-1] - avgY
        rotationVec[2i] = x[2i] - avgX
    end
    normalize!(rotationVec)
    d = dot(rotationVec,dx)
    dx .-= rotationVec .* d
    return nothing
end


# method 1
function genRandStrain!(dx::Vector{Float64},
                        x::Vector{Float64},
                        perturbMag::Number)
    # generate Strain 
    randn!(dx)
    normalize!(dx)
    dx .*= perturbMag

    # clean perturbation
    removeTranslation!(dx)
    removeRotation!(dx, x)
    return nothing
end

# method 2
function genRandStrain!(dx::Vector{Float64},
                        x::Vector{Float64},
                        perturbMag::Number,
                        nodes::Vector{Int})
    n = Int(length(x)/2)
    specialNodeMag = 2

    # generate Strain 
    randn!(dx)
    dx[ [2 .* nodes .- 1; 2 .* nodes] ] .*= specialNodeMag
    normalize!(dx)
    dx .*= perturbMag

    # clean perturbation
    removeTranslation!(dx)
    removeRotation!(dx, x)
    return nothing
end

# method 1
function genRandStrain(x::Vector{Float64},
                       perturbMag::Number
                       )::Vector{Float64}
    # add random perturbation to structure
    dx = similar(x)
    genRandStrain!(dx, x, perturbMag)
    return dx
end


# method 2
function genRandStrain(x::Vector{Float64},
                       perturbMag::Number,
                       nodes::Vector{Int}
                       )::Vector{Float64}
    # add random perturbation to structure
    dx = similar(x)
    genRandStrain!(dx, x, perturbMag, nodes)
    return dx
end




function applyRandStrain(x::Vector{Float64},
                         perturbMag::Number)::Vector{Float64}
    dx = genRandStrain(x, perturbMag)
    return  x .+ dx
end





######################
function addLig!(Q::Network, ligIndex::Int, siteIndex::Int)
    # ligIndex: 1=>solvent, 2=>right, 3=>wrong
    # siteIndex: 1=>active site, 2=>allosteric site 
    ligand = Q.ligs[ligIndex]
    site = Q.sites[siteIndex]
    @assert length(site) == length(ligand.bonds) 
    SM, RLM, EM = Q.pheno
    for (i, bond) in enumerate(ligand.bonds)
        c = site[i]
        SM[c] = SM'[c] = bond.k
        RLM[c] = RLM'[c] = bond.l
        EM[c] = EM'[c] = bond.e
    end
    return nothing
end


function removeLig!(Q::Network, ligIndex::Int, siteIndex::Int)

    ligand = Q.ligs[ligIndex]
    site = Q.sites[siteIndex]
    @assert length(site) == length(ligand.bonds) 
    SM, RLM, EM = Q.pheno
    for (i, bond) in enumerate(ligand.bonds)
        c = site[i]
        SM[c] = SM'[c] = 0.0
        RLM[c] = RLM'[c] = 0.0
        EM[c] = EM'[c] = 0.0
    end
    return nothing
end


function removeAllLig!(Q::Network)
    for i in eachindex(Q.pheno)
        Q.pheno[i] .*= Q.A
    end
    return nothing
end

######################################

function computeLengthMatrix(xy::Matrix{Float64},
                             A::Matrix{Float64})::Matrix{Float64}

    # Compute the length each bond in the network
    N = size(xy,1)
    lengthMat = zeros(size(A))
    @inbounds for i in 2:N
        for j in 1:i-1
            if A[i,j] == 0.
                nothing
            else
                deltax = xy[i,1] - xy[j,1]
                deltay = xy[i,2] - xy[j,2]
                lengthMat[i,j] = sqrt(deltax^2 + deltay^2)
            end
        end
    end
    lengthMat .= lengthMat .+ lengthMat'
    return lengthMat
end


function computeStress(xy::Matrix{Float64},
                       A::Matrix{Float64},
                       SM::Matrix{Float64},
                       RLM::Matrix{Float64})::Matrix{Float64}

    LM = computeLengthMatrix(xy, A)
    DM = SM .* (LM .- RLM)
    return DM
end



function computeStress(x::Vector{Float64},
                       A::Matrix{Float64},
                       SM::Matrix{Float64},
                       RLM::Matrix{Float64})::Matrix{Float64}
    numNodes = size(SM, 1)
    xy = x2xy(x)
    stressMatrix = computeStress(xy, A, SM, RLM)
    return stressMatrix   
end



function findStressedBond(stressMatrix::Array)
    bond = sortperm(stressMatrix[:])[1]
    bond = collect(Tuple(CartesianIndices(stressMatrix)[bond]))
    return bond
end


function findStressedNodes(stressMatrix::Matrix{Float64})::Vector{Int}

    bondCount = 15
    stressedBonds = sortperm(abs.(stressMatrix[:]))[end-bondCount:end]
    bonds_cart = CartesianIndices(stressMatrix)[stressedBonds]
    n = length(stressedBonds)
    nodes = zeros(2n)
    for i in 1:n
        nodes[i] = bonds_cart[i][1]
        nodes[i+n] = bonds_cart[i][2]
    end
    nodes = Int.(unique(nodes))
    return nodes
end


function constructParams(Q, assay)
    return Params(Q.unfoldedEnergy, Q.suppFitnessParams, assay)
end

function constructParamsArray(BigData, assay)
    paramsArray = Array{Params, ndims(BigData)}(undef, size(BigData))
    for i in eachindex(BigData)
        Q = BigData[i][1]
        paramsArray[i] = constructParams(Q, assay)
    end
    return paramsArray
end

function preligs2ligs(preligs::Vector{Vector{Vector{Float64}}})::Vector{Ligand}

    # This block converts the primitive ligand info
    # into structured ligand structs
    ligs = Ligand[]
    for i in 1:length(preligs) #Ligands
        listOfBondTrips = preligs[i]
        ligBondList = Bond[]
        for j in 1:length(listOfBondTrips)
            bondTrips = listOfBondTrips[j]
            push!(ligBondList, Bond(bondTrips...))
        end
        push!(ligs, Ligand(ligBondList))
    end
    return ligs
end

function ligs2preligs(ligs::Vector{Ligand})::Vector{Vector{Vector{Float64}}}
    # convert ligs to preligs        
    return [[ [bond.k, bond.l, bond.e] for bond in lig.bonds] for (i, lig) in enumerate(ligs)]
end


function makeSubDirs(dir::String, numSubDir::Int)
    # make subdirs within dir named 1,2,3 ...
    if dir[end] != '/'
        dir = dir * '/'
    end
    contents = readdir(dir)
    for i in 1:numSubDir
        if !isdir(dir*"$i")
            mkdir(dir*"$i")
        end
    end
end



function bondDataVector2Matrix(bondData::Vector{<:Number},
                               A::AbstractMatrix)
    @assert count(LowerTriangular(A) .!=0) == length(bondData)
    bondIndicies = findall(LowerTriangular(A).!=0)
    bondDataMatrix = zeros(size(A))
    bondDataMatrix[bondIndicies] .= bondData
    symmetric!(bondDataMatrix)    
    return bondDataMatrix
end

function convertSymMat2Vec(M::Matrix)
    
    # store the lower triangular non-zero entries
    # of a symetric square matrix in a vector.
    
    @assert issymmetric(M)
    return sparse(LowerTriangular(M)).nzval
end

function convertVec2SymMat(A::Matrix, v::Vector)
    # Creat matrix M and fill it with the
    # nonzero positions of A with values of v.
    
    @assert issymmetric(A)
    # make sure the length of v matches the number of non zero
    # entries in the lower triangle of A.
    @assert sum( LowerTriangular(A) .!= 0 ) == length(v) 
    
    M = zeros(eltype(v), size(A))
    M[findall(!iszero, LowerTriangular(A))] .= v
    symmetric!(M)
    return M
end
