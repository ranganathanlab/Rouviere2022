module CoreFuncs

using LinearAlgebra
using Random
using Statistics
using Distributed
using Printf
using SparseArrays

push!(LOAD_PATH,"./")

#import General: eye, acceptMC, symmetric, symmetric!, computeProgress, intol


include("constants.jl")
include("CoreFuncs/types.jl")
include("CoreFuncs/model.jl")
include("CoreFuncs/relax.jl")
include("CoreFuncs/mutate.jl")
include("CoreFuncs/fitnessfuncs.jl")
include("CoreFuncs/tools.jl")
include("CoreFuncs/groundstate.jl")
include("CoreFuncs/buildnetwork.jl")
include("CoreFuncs/evolve.jl")

export

    # from types.jl
    Bond,
    Ligand,
    Network,
    Params,
    LilNetwork,
    LiteNetwork,

    # from buildnetwork.jl
    buildNetwork,
    buildNetworkFromSettings,
    buildxy0,
    convertNetwork2LiteNetwork,
    convertLiteNetwork2Network,

    # from evolve.jl
    evolve,
#    evolve2,
    evolveLite,
    evolveFlux,
    runMany,
    runManyFlux,
    evolveFamily,
    formatFamily,
    makeEvoTrace,
    setLigs!,
    evolvePop,
    pickWinners,
    pickLosers,
    evolvePopWrapper,

    # from fitnessfuncs.jl

    stability,
    instability,
    binding,
    specificity,
    specificity2,
    specificity3,
    stabilityBinding,
    stabilitySpecificity,
    allostery,
    allostery2,
    computeFitness,
    setLigDex,
    computeProbs,
    PAllostery3,
    conc2pot,
    custom,

    # from groundstate.jl
    findGroundState,
    findGroundStates,
    findGroundStates!,
    findGroundStates2,
    swapLigs,
    findStates,
    buildRegistry,
    makeRegistryUnique,
    removeDuplicates,

    # from model.jl 
    computeEnergy,
    computeForce,
    computeGradient!,
    computeHessian,
    computeHessian!,

    # from mutate.jl
    mutateBond!,
    mutateNode!,
    randomMutateNode!,
    mutateNodeCouple!,
    mutate!,

    # from relax.jl
    relaxSprings,
    minimizeFIRE,
    implicitEuler!,
    velocityVerlet!,
    relaxSpringsGD!,
    relaxSpringsGD,

    # from tools.jl
    xy2x,
    x2xy,
    computeDistanceMatrix,
    sameConf,
    findLowEIndices,
    genRandStrain!,
    genRandStrain,
    applyRandStrain,
    center!,
    measureAngle,
    rotate!,
    orient!,
    orientNetwork!,
    addLig!,
    removeLig!,
    removeAllLig!,
    computeStress,
    computeConfDiff,
    constructParams,
    constructParamsArray,
    preligs2ligs,
    ligs2preligs,
    makeSubDirs,
    bondDataVector2Matrix,
    convertSymMat2Vec,
    convertVec2SymMat,
    intol
end # end of module

