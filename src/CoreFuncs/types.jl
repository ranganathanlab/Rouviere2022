# define the main type of the program


mutable struct Bond
    k::Float64 # spring constant
    l::Float64 # restlength
    e::Float64 # energy
end

mutable struct Ligand
    bonds::Vector{Bond} # bond
end


mutable struct Network 
    A::Matrix{Float64} 
    NRLM::Matrix{Float64}
    geno::Vector{Int64} 
    pheno::Vector{Matrix{Float64}} 
    typeTable::Array{Float64,3} 
    energies::Matrix{Float64} 
    structs::Matrix{Matrix{Float64}}
    ligs::Vector{Ligand}
    sites::Vector{Vector{CartesianIndex{2}}}
    unfoldedEnergy::Float64
    suppFitnessParams::Tuple
end



struct Params
    unfoldedEnergy::Float64
    suppFitnessParams::Tuple
    assay::String
end


struct LilNetwork
    # stores only the things that are different across networks of an ensemble.
    geno::Vector{Int64}
    energies::Matrix{Float64}
    structs::Matrix{Matrix{Float64}}
end

struct LiteNetwork
    A::BitMatrix
    NRLM::Vector{Float64}
    geno::Vector{Int64} 
    typeTable::Array{Float64,3} 
    energies::Vector{Float64} 
    structs::Vector{Matrix{Float64}}
    preligs::Vector{Vector{Vector{Float64}}}
    sites::Vector{Vector{CartesianIndex{2}}}
    unfoldedEnergy::Float64
    suppFitnessParams::Tuple
    ligDex::Vector{CartesianIndex{2}}
end


