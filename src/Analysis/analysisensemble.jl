#######################################################
########## Ensemble Analysis Functions #########
#######################################################


function ensembleComputeGSEMutSens(ensemble,
                                   ligDex::Vector{CartesianIndex{2}},
                                   mutType::String;
                                   saveStructs::Bool=false,
                                   GS_settings::Tuple{Int, Int}=(10,100),
                                   bondStrain::Number=0.2)
    # Computes deep mutational scan
    s = size(ensemble)
    l = ndims(ensemble)
    outList = pmap(x -> computeGSEMutSens(x, ligDex, mutType; GS_settings, saveStructs, bondStrain), ensemble)

    if saveStructs
        GSE_WT_ens = Array{typeof(outList[1][1]),l}(undef,s)
        GSE_DMS_ens = Array{typeof(outList[1][2]),l}(undef,s)
        structs_DMS_ens = Array{typeof(outList[1][3]),l}(undef,s)
        for i in eachindex(outList)
            GSE_WT_ens[i] = outList[i][1]
            GSE_DMS_ens[i] = outList[i][2]
            structs_DMS_ens[i] = outList[i][3]
        end
        return GSE_WT_ens, GSE_DMS_ens, structs_DMS_ens
    else
        GSE_WT_ens = Array{typeof(outList[1][1]),l}(undef,s)
        GSE_DMS_ens = Array{typeof(outList[1][2]),l}(undef,s)
        for i in eachindex(outList)
            GSE_WT_ens[i] = outList[i][1]
            GSE_DMS_ens[i] = outList[i][2]
        end
        return GSE_WT_ens, GSE_DMS_ens
    end
end



function ensembleComputeFitSens(GSE_WT_ens, GSE_DMS_ens, params::Params)
    return map((x,y) -> computeFitSens(x, y, params), GSE_WT_ens, GSE_DMS_ens)
end


function ensembleComputeFitness(GSE_WT_ens, params::Params)
    return map(x -> computeFitness(x, params), GSE_WT_ens)
end


function ensembleComputeSurfaceAllostery(GSE_WT_ens, GSE_DAS_ens, params::Params, sign::String)
    return map( (x,y) -> computeSurfaceAllostery(x, y, params, sign), GSE_WT_ens, GSE_DAS_ens)
end


function ensembleConfChangeTest(ensemble, ligDex)
    return map(x -> confChangeTest(x, ligDex), ensemble)
end

function ensembleConfChangeTest(ensemble, lig1::Ligand, lig2::Ligand, ligDex)
    return map(x -> confChangeTest(x, lig1, lig2, ligDex), ensemble)
end

function ensembleLocateConfChange(ensemble, ligDex)
    return map(x -> locateConfChange(x, ligDex), ensemble)
end


function ensembleConfChangeTest2(ensemble, ligDex)
    return map(x -> confChangeTest2(x, ligDex), ensemble)
end

function ensembleComputeEnergyLigSlice(ensemble, ligLengths::AbstractVector{<:Number})
    return pmap( x -> computeEnergyLigSlice( x, ligLengths), ensemble )   
end


function ensembleComputeLigScape(ensemble, ligLengths; GS_settings=(10,100))
    out = pmap( x -> computeLigScape(x, ligLengths, GS_settings=GS_settings), ensemble )
    energies_ens = Array{typeof(out[1][1]),ndims(ensemble)}(undef, size(ensemble))
    structs_ens = Array{typeof(out[1][2]), ndims(ensemble)}(undef, size(ensemble))
    for i in eachindex(out)
        energies, structs = out[i]
        energies_ens[i] = energies
        structs_ens[i] = structs
    end
    return energies_ens, structs_ens
end

function ensembleAverageBondStrain(ensemble::Vector{Network}, ligDex::Vector{CartesianIndex{2}})
    @assert length(ligDex) == 2
    cummBondStrains = zeros(size(ensemble[1].A))
    for i in eachindex(ensemble)
        Q = ensemble[i]
        cummBondStrains .+= abs.(computeBondStrains(Q.structs[ligDex[1]], Q.structs[ligDex[2]], Q.A))
    end
    return cummBondStrains ./ length(ensemble)
end 

function ensembleComputeOverlaps(ensemble::Array{Network})
    return computeOverlaps.(ensemble)
end


function ensembleComputeFracOfSites(df_ens, sign, thresh)
    return map(x -> computeFracOfSites(x, sign, thresh), df_ens)
end









