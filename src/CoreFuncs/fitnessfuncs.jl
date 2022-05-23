
stability(E_s) = -E_s
specificity(E_s, E_r, E_w) = min(E_s-E_r, E_w-E_s)
binding(E_s, E_r) = E_s - E_r
doubleBinding(E_s, E_r, E_w) = min(E_s-E_r, E_s-E_w)
allostery(E00, E10, E01, E11) = (E10-E00) - (E11-E01) 
negativeAllostery(E00, E10, E01, E11) = -allostery(E00, E10, E01, E11)

function customSpecificity(E_s, E_r, E_w, d_r, d_w, α)
    # a function that returns how closely the three
    # ground state energies, E_s, E_r, E_w match the 
    # target energies differences  d_r, d_w.
    # α give the tolerance of mathcing.
    

    dE_r = E_r - E_s
    dE_w = E_w - E_s
    s = abs(dE_r - d_r) + abs(dE_w - d_w)
    return -s * α
end

function setLigDex(assay::String, ligSpaceSize=(5,5))
    # returns the indices of the ligands senarors requried by
    # by a assay from the argument.

    C = CartesianIndices(ligSpaceSize)
    
    if assay in ["Stability", "solvent"]
        c = [C[1,1]]
    elseif assay in ["Binding", "Big2", "sr"]
        c=[C[1,1], C[2,1]]
    elseif assay in ["Specificity", "DoubleBinding", "CustomSpecificity"]
        c=[C[1,1],C[2,1],C[3,1]]
    elseif assay in ["Allostery","NegativeAllostery","Big4"]
        c=[C[1,1],C[2,1],C[1,2],C[2,2]]
    elseif assay=="Big2Allo"
        c=[C[1,1],C[2,2]]
    elseif assay=="Big5"
        c=[C[1,1],C[2,1],C[3,1],C[1,2],C[2,2]]
    elseif assay=="All"
        c=C[:]
    elseif assay=="sw"
        c=[C[1,1], C[3,1]]
    elseif assay=="rw"; c=[2,3]
        c=[C[2,1], C[3,1]]
    elseif assay in ["SurfaceAllostery", "NegativeSurfaceAllostery"]
        c=[C[1,4], C[2,4], C[1,5], C[2,5]]
    else
        error("You spelled wrong idiot!")
    end
    return c
end



function fitness(energies::Matrix{Float64},
                 unfoldedEnergy::Float64,
                 suppFitnessParams::Tuple,
                 assay::String)::Float64

    if assay == "Stability"
        f = stability(energies[setLigDex(assay)]...)
    elseif assay == "Binding"
        f = binding(energies[setLigDex(assay)]...)
    elseif assay == "Specificity"
        f = specificity(energies[setLigDex(assay)]...)
    elseif assay == "DoubleBinding"
        f = doubleBinding(energies[setLigDex(assay)]...)
    elseif assay == "Allostery"
        f = allostery(energies[setLigDex(assay)]...)
    elseif assay == "NegativeAllostery"
        f = negativeAllostery(energies[setLigDex(assay)]...)
    elseif assay == "SurfaceAllostery"
        f = allostery(energies[setLigDex(assay)]...)
    elseif assay == "NegativeSurfaceAllostery"
        f = negativeAllostery(energies[setLigDex(assay)]...)
    elseif assay == "CustomSpecificity"
        f = customSpecificity(energies[setLigDex(assay)]..., suppFitnessParams...)
    else
        error("You spelled wrong idiot!")
    end

    # Add Stability pressure
    E_s = energies[setLigDex("solvent")...]
    E_s > unfoldedEnergy && (f = -10 - E_s)
    return f
end


function computeFitness(energies::Matrix{Float64},
                        unfoldedEnergy::Float64,
                        suppFitnessParams::Tuple,
                        assay::String)::Float64

    ligDex = setLigDex(assay) 
    return fitness(energies, unfoldedEnergy, suppFitnessParams, assay)
end

function computeFitness(Q::Network, assay::String)::Float64
    return computeFitness(Q.energies, Q.unfoldedEnergy, Q.suppFitnessParams, assay)
end

function computeFitness(energies::Matrix{Float64}, params::Params)::Float64
    return computeFitness(energies, params.unfoldedEnergy, params.suppFitnessParams, params.assay)
end




