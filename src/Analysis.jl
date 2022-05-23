
module Analysis

using Statistics
using LinearAlgebra
using Distributed
using Printf
push!(LOAD_PATH, ".")
#import General: linspace, eye, l2a, k2ij, vec2UpperLowerTriMat!, intol
using CoreFuncs


include("constants.jl")

# analysis functions
include("Analysis/mutscan.jl")
include("Analysis/analysis.jl")
include("Analysis/evolvability.jl")
include("Analysis/analysisensemble.jl")
include("Analysis/epistasis.jl")
include("Analysis/families.jl")
export

    # from mutScan
    computeGSEMutSens,
    computeGSEMutSens,
    computeFitSens,
    computeSingleMutsFit,
    doMutScan!,
    testAllMutantsWithAllStructures!,
    removeRedundantStructures,
    getMutDetailsList,
    computeSurfaceAllostery,
    findAllostericHotSpot,
    getSurfaceBonds,

    # from analysis.jl
    extractEnergyTrace,    
    computeStatehoodAllo,
    confChangeTest,
    confChangeTest2,
    locateConfChange,
    computeBondStrains,
    computeEnergyLigSlice,
    computeLigScape,
    labelLigScape,
    outlineLigScape,
    computeLeastStepMatrix,
    computeStepsFromAct,
    getUniqueSingleMutStructs,
    doubleGaussian,
    computeBimodality,
    getDeformation,           
    computeOverlaps,
    computeOverlapsStiffnesses,
    measureInterFrac, 
    countSoftBonds,
    computeFracOfSites,
    computeCCSize,
    computeStiffnessOfAlloMode1,
    computeStiffnessOfAlloMode2,
    measureBarrier, 
    computeEnergyAlongCoor,
    computeRuggedness, 

    # from analysisensemble.jl
    ensembleComputeGSEMutSens,
    ensembleConfChangeTest,
    ensembleConfChangeTest2,
    ensembleLocateConfChange,
    ensembleComputeSurfaceAllostery,
    ensembleComputeFitSens,
    ensembleComputeFitness,
    ensembleComputeEnergyLigSlice,
    ensembleComputeLigScape, 
    ensembleAverageBondStrain,
    ensembleComputeGSEBondSens,
    ensembleComputeOverlaps,
    ensembleComputeFracOfSites,

    # from epistasis.jl
    computeDoubleMutGSE,
    doubleMutGSE2Fit,
    computeEpistasis,
    computeEpistasisArray,
    makeEpiMat,   
    fastTrackEpiMat,
    
    # from evolvability.jl
    computeEvoTraces,
    ensembleComputeEvoTraces,

    # from families.jl
    makeMSA,
    computeSeqSim,
    computeSeqSimMat,
    dec2bin,
    computeFreqs,
    computeConservation,
    relEnt,
    computeWeights,
    computeCorr,
    compressCorr

end
