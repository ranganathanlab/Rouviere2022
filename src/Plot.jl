
module Plot

using Statistics
using LinearAlgebra
using Distributed
using Printf
using PyCall
using PyPlot
using LaTeXStrings
using Distributions
push!(LOAD_PATH, ".")

#import General: linspace, eye, l2a, k2ij, vec2UpperLowerTriMat!
using CoreFuncs
using Analysis

include("constants.jl")

# plotting functions
include("Plot/plottingtools.jl")
include("Plot/ploto.jl")
include("Plot/singlefigures.jl")
include("Plot/ensemblefigures.jl")
include("Plot/makesinglefigs.jl")


export
    
    # from ploto.jl
    plotNetwork,
    plotNetwork!,
    plotComplex!,
    plotComplex,
    plotBondStrains!,
    plotCorrelations,
    plotScatter!,
    wilsonCIHelper,
    waldCIHelper,

    # from singlefigures.jl
    genEnergiesFig,
    genEnergiesAlloFig,
    genStructsFig,
    genStructsAlloFig,
    gen3SFig,
    genLowEnergyTableFig,
    genMutSensMapFig,
    genMutSensMapAbsFig,
    genMutSensMapRawFig,
    genDMSFig,
    genMutSensDistFig,
    genMutSensRDFFig,
    genStrainMapFig,
    genStrainAlloMapFig,
    genLigEnergySliceFig,
    genLigScapeFig,
    genEnergySpectrumFig,
    genEpiMatrixFig,
    genEpiHistFig,
    genEpiMatSpectrumFig,
    genEvolFig, 
    genStabFitDMSFig,
    genStabBindFig,
    genStabBindScatterHistFig,
    genEnergyTickFig,
    genSeqSimFig,
    genColoredStructFig, 
    genConsMapFig,
    genPosMatrixFig,
    genFitnessesFig,
    genBondSensMapFig,
    genStrainVsSoftModeFig,
    genRegistryFig,
    #genLigScapeAlloLabelStatesFig,
    genConfChangeDMSFig,
    genEnergyVsConfCoorFig,
    genEnergyVsConfCoorFig2,
    genSurfaceAllosteryFig,


    # from ensemblefigures.jl
    genConfChangeBarFig,
    genLocConfChangeBarFig,
    genAlloBarFig,
    genAlloFitFig,
    genAlloCCLFig,
    genCCDFitFig,
    genStabFitFig,
    gen_q1_vs_q234Fig,

    # from makessinglefigs.jl
    makeAll_energies,
    makeAll_energyVsConfCoor,
    makeAll_structs,
    makeAll_structsAllo,
    makeAll_strain,
    makeAll_strainAllo,
    makeAll_mutSensMap,
    makeAll_DMS,
    makeAll_mutSensDist,
    makeAll_mutSensRDF,
    makeAll_StabBind,
    makeAll_mutSens,
    makeAll_ligScape, 
    makeAll_strainVsSoftMode,

    # from plottingtools.jl
    formatFig!,
    jitter!,
    makeAllNetworkAxSame!,
    makeColorbar!
end
