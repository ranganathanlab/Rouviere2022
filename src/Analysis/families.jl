
function makeMSA(listOseqs::Vector{Vector{Int}})
    numSeqs = length(listOseqs)
    seqLen = length(listOseqs[1])

    msa = zeros(Int, seqLen, numSeqs)
    for i in 1:numSeqs
        msa[:,i] .= listOseqs[i]
    end
    return msa
end


function dec2bin(msa::Matrix{Int}; libSize::Int=5)
    # convert decimal msa to binary msa
    seqLen, numSeqs = size(msa)
    X = zeros(seqLen*libSize, numSeqs)
    for j in 1:numSeqs
        seq = msa[:,j]
        for i in 1:seqLen
            index = (i-1) * libSize + seq[i]
            X[index,j] = one(Float64)
        end
    end
    return X
end


function computeSeqSim(seq1::AbstractVector{<:Number},
                       seq2::AbstractVector{<:Number})
    N = length(seq1)
    t = 0::Int
    @inbounds @simd for i in 1:N
        t += ifelse(seq1[i]==seq2[i],1,0)
    end
    return t / N
end

function  computeSeqSimMat(msa::Matrix{Int})
    # msa has sequences as cols
    seqLen, numSeqs = size(msa)
    seqSimMat = zeros(numSeqs, numSeqs)
    for j in 1:numSeqs, i in 1:numSeqs
        @views seqSimMat[i,j] = computeSeqSim(msa[:,j], msa[:,i])
    end
    return seqSimMat
end

function computeSeqSimMat(X::Matrix{Float64};
                          libSize::Int=5)::Matrix{Float64}
    a = libSize / size(X,1)
    return (X' * X ) .* a
end


# not tested
function computeFreqs(X::Matrix{Float64};
                      libSize::Int=5,
                      f0::Vector{Float64}=ones(libSize) ./ libSize,
                      λ::Number=0.001)
#    @assert length(f0) == libSize
    numSeqs=size(X,2)
    seqLen = Int(size(X,1) / libSize)
    w = ones(numSeqs) ./ numSeqs # seq weights
    f1 = X * w
    f2 = (X .* w') * X'
    f2_bkg =  repeat(f0 * f0', seqLen, seqLen)
    for i in 1:seqLen
        inds = libSize*(i-1)+1 : libSize*i
        f2_bkg[inds, inds] = diagm(f0)
    end
    f1_reg = (1 - λ) .* f1 .+ λ .* repeat(f0, seqLen)
    f2_reg = (1 - λ) .* f2 .+ λ .* f2_bkg
  #  f0_reg = mean(reshape(f1_reg, seqLen, libSize), dims=1)
   # return f1_reg, f2_reg  #, f0_reg
    return f1_reg, f2_reg
end

# not tested
function relEnt(P::AbstractVector{T}, Q::AbstractVector{T}) where T<:Number
    # compute relative entropy of p compared to q
    @assert length(P) == length(Q)
    P .= P ./ sum(P)
    Q .= Q ./ sum(Q)
    D = zero(Float64)
    @inbounds for i in eachindex(P)
        p = P[i]
        q = Q[i]
        r = p/q
        if q != zero(eltype(Q)) && p != zero(eltype(Q))
            D += p*log(p/q) 
        end
    end
    return D
end

# not tested
function computeConservation(f1::Vector{Float64};
                             libSize::Int=5,
                             f0::Vector{Float64}=ones(libSize) ./ libSize)
    # compute positional conservation by relative entropy
    @assert length(f0) == libSize
    seqLen = Int(length(f1)/libSize)
    f1_  = reshape(f1, :, seqLen)'
    D = zeros(seqLen)
    for i in eachindex(D)
        D[i] = relEnt(f1_[i,:],f0)
    end
    return D
end


function computeCorr(f1::Vector{Float64},
                     f2::Matrix{Float64})
    return f2 .- (f1 * f1')
end


function computeWeights(f1::Vector{Float64};
                        libSize::Int=5,
                        f0::Vector{Float64}=ones(libSize) ./ libSize)

    seqLen = Int(length(f1) / libSize) 
    ϕ = zeros( size(f1) )
    for i in 1:seqLen, j in 1:libSize 
        q = f0[j]
        k = libSize*(i-1) + j
        f = f1[k]
        ϕ[k] = log( (f*(1-q)) / ((1-f)*q) )
    end
    return ϕ
end

function compressCorr(C::Matrix{Float64}; libSize::Int=5)

    seqLen = Int(size(C,1)/libSize)
    Cij = zeros(seqLen, seqLen)
    for i in 1:seqLen, j in 1:seqLen
        inds_i = libSize*(i-1) + 1 : libSize*i
        inds_j = libSize*(j-1) + 1 : libSize*j
        Cij[i,j] = norm(C[inds_i, inds_j][:])
    end
    return Cij

end


function mutInfMat(f2::Matrix, f1::Vector; libSize=5)
    
    seqLen = Int(length(f1) / libSize)
    M = zeros(seqLen,seqLen)
    f1f1 = f1 * f1'

    for i in 1:seqLen, j in 1:seqLen
        p2 = f2[ 5*(i-1)+1:5*i, 5*(j-1)+1:5*j ][:]
        p1 = f1f1[ 5*(i-1)+1:5*i, 5*(j-1)+1:5*j ][:]
        p1 .= p1 ./ sum(p1)
        p2 .= p2 ./ sum(p2)
        M[i,j] = relEnt(p2,p1)
    end
    return M
end

export mutInfMat




















