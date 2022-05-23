

function genEnergiesFig(fitnesses::Vector{Float64}, energies::Matrix{Float64})
    # Plots the energy of sol, right and wrong  ligand bound cases
    # returns figure object
    numMuts = length(fitnesses)
    energies = energies'
    fig, ax = subplots(2,1, figsize=(10,6))
    formatFig!(fig)
    ax[1].plot(1:numMuts, energies[:,1], c=clr("blue"))
    ax[1].plot(1:numMuts, energies[:,2], c=clr("green"))
    ax[1].plot(1:numMuts, energies[:,3], c=clr("red"))
    ax[1].set_xlabel("mutation number")
    ax[1].set_ylabel("energy")
    ax[1].legend(["solvent","right","wrong"])
    
    ax[2].plot(1:numMuts, fitnesses, c=clr("blue"))
    ax[2].plot(1:numMuts, zeros(numMuts), c=clr("black"), lw=1)
    ax[2].set_xlabel("mutation number")
    ax[2].set_ylabel("fitness")

    if numMuts < 25
        ax[1].set_xticks(1:numMuts)
        ax[2].set_xticks(1:numMuts)
    end

    fig.subplots_adjust(hspace=0.5)
    return fig
end


function genEnergiesAlloFig(fitnesses::Vector{Float64}, energies::Matrix{Float64})
    numMuts = length(fitnesses) 
    energies = energies'

    fig, ax = subplots(2,1, figsize=(10,6))
    formatFig!(fig)
    ax[1].plot(1:numMuts, energies[:,2] .- energies[:,1], c=clr("blue"))
    ax[1].plot(1:numMuts, energies[:,4] .- energies[:,3], c=clr("green"))

    ax[1].set_xlabel("mutation number")
    ax[1].set_ylabel("energy")
    ax[1].legend(["\$ ΔE_{na}\$","\$ ΔE_a\$"])

    ax[2].plot(1:numMuts, fitnesses, c=clr("blue"))
    ax[2].plot(1:numMuts, zeros(numMuts), c=clr("black"), lw=1)
    ax[2].set_xlabel("mutation number")
    ax[2].set_ylabel("fitness")
    
    if numMuts < 25
        ax[1].set_xticks(1:numMuts)
        ax[2].set_xticks(1:numMuts)
    end

    fig.subplots_adjust(hspace=0.5)
    return fig
end

function genEnergyTickFig(Q::Network)
    f!(x, y, color) = ax.scatter([x], [y], 1000, marker="_", c=color )
    fig, ax = subplots(figsize=(4,6))
    formatFig!(fig)
    fig.subplots_adjust(left=0.25)
    blue = [0.2 .2 .8]
    red = [.8 .2 .2]
    green = [0.2 0.6 0.2]
    xlim = [-0.3, 2.3]
    f!(0, Q.energies[1], blue )
    f!(0, Q.energies[2], green )
    f!(0, Q.energies[3], red)
    f!(1, Q.energies[4], blue )
    f!(1, Q.energies[5], green)
    f!(1, Q.energies[6], red)
    f!(2, Q.energies[7], blue )
    f!(2, Q.energies[8], green)
    f!(2, Q.energies[9], red)
    ax.plot(xlim, Q.unfoldedEnergy*ones(size(xlim)), c="k", ls="--")
    ax.set_xlim(xlim)
    ax.set_ylim([0, 1.05*1*max(Q.unfoldedEnergy, maximum(Q.energies))])
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_ylabel("Energy", fontsize=16)
    return fig, "energies"
end

function genStructsFig(Q::Network)

    numNodes = size(Q.A, 1)
    SM, RLM, EM = Q.pheno
    structs = Q.structs[[1,2,3]] # get only good ones
    numStructs = length(structs)
    orientNetwork!(Q)
    fig, ax = subplots(1,numStructs, figsize=(15,5))
    formatFig!(fig); fs = 20
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)

    colors=["b","g","r"]
    titles=["Solvent","Right","Wrong"]
    for i in 1:numStructs
        plotNetwork!(ax[i], structs[i], SM, bondColor=[.5,.5,.5])
        ax[i].axis("off")
        ax[i].axis("equal")

        # plot active site ligand
        actSite = Q.sites[1]
        for (j,loc) in enumerate(actSite)
            ax[i].plot( [structs[i][loc[1],1], structs[i][loc[2],1]],
                       [structs[i][loc[1],2], structs[i][loc[2],2]],
                       color=colors[i], linewidth=3, zorder=1 )
        end

        # plot active site ligand
        alloSite = Q.sites[2]
        for (j,loc) in enumerate(alloSite)
            ax[i].plot( [structs[i][loc[1],1], structs[i][loc[2],1]],
                       [structs[i][loc[1],2], structs[i][loc[2],2]],
                       color=colors[1], linewidth=3, zorder=1 )
        end

        ax[i].set_title(titles[i], y=1,fontsize=fs)
    end
    makeAllNetworkAxSame!(ax)
    return fig, "structs"
end


function genStructsAlloFig(Q::Network; bondWidth=3, nodeSize=30, assay="Allostery")

    numNodes = size(Q.A, 1)
    numLigs = length(Q.ligs)
    SM, RLM, EM = Q.pheno

    ligDex = setLigDex(assay)
    structs = Q.structs[ligDex]
    orientNetwork!(Q)
    fig, ax = subplots(2,2,figsize=(10,8))
    formatFig!(fig); fs = 16
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)
    
    #colors=["b", "g"]
    colors=[[0,0,1], [0.1, .65 ,0.1]]
    actSiteColors = [1,2,1,2]
    alloSiteColors = [1,1,2,2]
    titles = ["00","10","01","11"]
    for i in eachindex(structs)
        xy = structs[i]
        plotNetwork!(ax[i], xy, SM, bondColor=[.5,.5,.5], bondWidth=bondWidth, nodeSize=nodeSize)
        ax[i].axis("off")
        ax[i].axis("equal")
        ax[i].set_title(titles[i], fontsize=20)
        # plot active site ligand
        actSite = Q.sites[1]
        for (j,loc) in enumerate(actSite)
            ax[i].plot( [xy[loc[1],1], xy[loc[2],1]],
                       [xy[loc[1],2], xy[loc[2],2]],
                       color=colors[actSiteColors[i]], linewidth=bondWidth+2, zorder=1 )
        end

        # plot active site ligand
        alloSite = Q.sites[2]
        for (j,loc) in enumerate(alloSite)
            ax[i].plot( [xy[loc[1],1], xy[loc[2],1]],
                       [xy[loc[1],2], xy[loc[2],2]],
                       color=colors[alloSiteColors[i]], linewidth=bondWidth+2, zorder=1 )
        end

    end
    makeAllNetworkAxSame!(ax)
    return fig, "structs"
end



function genSurfaceAllosteryFig(Q::Network, df::Matrix{Float64})

    titles = ["Negative Bond Strain", "Positive Bond Strain"]
    cmap = "Reds"
    numShades = 100
    cmp = get_cmap(cmap,numShades)
    numNodes = size(Q.A, 1)
    surfaceBonds = getSurfaceBonds()
    @assert length(surfaceBonds)  == size(df,1)
    SM, RLM, EM = Q.pheno
    orientNetwork!(Q)
    xy = Q.structs[1] # get only good ones

    fig, axx = subplots(1,2, figsize=(10,5))
    formatFig!(fig); fs = 20
    #fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)
    
    for i in 1:length(axx)
        ax = axx[i]
        colors=["b","g","r"]
        plotNetwork!(ax, xy, SM, bondColor=[.5,.5,.5])
        ax.axis("off")
        ax.axis("equal")

        # plot active site ligand
        actSite = Q.sites[1]
        for (j,loc) in enumerate(actSite)
            ax.plot( [xy[loc[1],1], xy[loc[2],1]],
                       [xy[loc[1],2], xy[loc[2],2]],
                       color="b", linewidth=3, zorder=1 )
        end
        
        # plot active site ligand
        alloSite = Q.sites[2]
        for (j,loc) in enumerate(alloSite)
            ax.plot( [xy[loc[1],1], xy[loc[2],1]],
                       [xy[loc[1],2], xy[loc[2],2]],
                       color="b", linewidth=3, zorder=1 )
        end
        
        df1 = abs.(df[:,i])
        df1 = Int.(round.(df1.* ( numShades/ maximum(abs.(df)))))
        c = cmp.(df1)

        for (i, bond) in enumerate(surfaceBonds)
            ax.plot( [xy[bond[1],1], xy[bond[2],1]],
                       [xy[bond[1],2], xy[bond[2],2]],
                       color=c[i], linewidth=3, zorder=1 )
        end
        ax.set_title(titles[i])
    end 

    return fig, "SSA"
end



function genStrainMapFig(Q::Network; colorTop::Number=0.2)

    fig, ax = subplots(1,3, figsize=(15,5))
    formatFig!(fig)
    fig.subplots_adjust(left = 0.1,bottom=0.1, top=0.9, right=0.85)
    fs = 16
    pairs = [[1,2],[1,3],[2,3]]
    labes = ["Solvent", "Right", "Wrong"]
    orientNetwork!(Q)
    removeAllLig!(Q)
    SM, RLM, EM = Q.pheno
    for i in 1:size(ax,1)
        first = pairs[i][1]
        second = pairs[i][2]
        strainMatrix = computeBondStrains(Q.structs[first], Q.structs[second], SM)
        plotColoredBonds!(ax[i], Q.structs[1], SM, strainMatrix, colorTop)
        ax[i].text(-0.05, 0.5, "$(labes[first])-$(labes[second])",
            transform=ax[i].transAxes, rotation=90, fontsize=16,
            horizontalalignment="center", verticalalignment="center")
    end
    makeAllNetworkAxSame!(ax)
    makeColorbar!(fig, cm="Greens", low=0, high=colorTop, x=.9, y=.25, wid=.02, len=.5)
    fig.text(.94, .50, "Strain", rotation=270, fontsize=16,
         horizontalalignment="center", verticalalignment="center")
    return fig, "strainMap"
end



function genStrainAlloMapFig(Q::Network; colorTop::Number=0.2, assay="Allostery")

    fig, ax = subplots(1,3,figsize=(15,5))
    formatFig!(fig)
    fig.subplots_adjust(left = 0.1,bottom=0.1, top=0.9, right=0.85)
    fs = 16
    titles = ["Before", "After"]
    
    ligDex = setLigDex(assay)
    




    pairs = [[ligDex[1] ligDex[2]], [ligDex[1],ligDex[3]], [ligDex[1],ligDex[4]]]
    labes = ["00-10", "00-01", "00-11"]

    orientNetwork!(Q)
    SM, RLM, EM = Q.pheno
    for i in 1:length(ax)
        first = pairs[i][1]
        second = pairs[i][2]
        strainMatrix = computeBondStrains(Q.structs[first], Q.structs[second], SM)
        plotColoredBonds!(ax[i], Q.structs[ligDex[1]], SM, strainMatrix, colorTop)
        ax[i].text(0.5, 1.1, labes[i],
               transform=ax[i].transAxes, fontsize=16,
               horizontalalignment="center", verticalalignment="center")
    end
    makeAllNetworkAxSame!(ax)
    makeColorbar!(fig, cm="Greens", low=0, high=colorTop, x=.9, y=.25, wid=.02, len=.5)
    fig.text(.94, .50, "Strain", rotation=270, fontsize=16,
        horizontalalignment="center", verticalalignment="center")
    return fig, "strainMap"
end



function genMutSensMapFig(Q::Network,
                          df::Matrix{Float64},
                          assay::String,
                          cmap::String,
                          vmax::Number,
                          vmin::Number)

    # mutational sensitivity map
    fig, ax = subplots(figsize=(5,4))
    formatFig!(fig)
    fig.subplots_adjust(left = 0.1,bottom=0.1, top=0.9, right=0.90)
    fs = 16
    orientNetwork!(Q)
    df = vec(mean(df,dims=2))
    plotColoredStruct!(ax, Q.structs[1], Q.A, df, cmap=cmap, vmax=vmax, vmin=vmin)

    ax.axis("equal")
    ax.axis("off")
    makeColorbar!(fig, cm=cmap, low=vmin, high=vmax, x=.92, y=.20, wid=.04, len=.6)
    
    fig.text(.99, .5, "Δ $assay", rotation=270, fontsize=16,
             horizontalalignment="center", verticalalignment="center")
    return fig, "mutSensMap_" * assay
end



function genMutSensMapAbsFig(Q::Network,
                             df::Matrix{Float64},
                             assay::String;
                             cmap="Blues",
                             vmax=maximum(mean(abs.(df),dims=2)),
                             vmin=0)
    # plot the absolute value of mutational sensitivity on the strucuture.
    fig, name = genMutSensMapFig(Q, abs.(df), assay, cmap, vmax, vmin)
    return fig, name
end

function genMutSensMapRawFig(Q::Network,
                             df::Matrix{Float64},
                             assay::String;
                             cmap="bwr",
                             vmax=maximum(abs.(df)) /2,
                             vmin=-maximum(abs.(df)) /2 )
    # plot the raw value of mutational sensitivity on the strucuture.
    fig, name = genMutSensMapFig(Q, df, assay, cmap, vmax, vmin)
    return fig, name
end



function genBondSensMapFig(Q::Network,
                           df::Matrix{Float64},
                           assay::String;
                           cmap="Blues",
                           vmax=maximum(abs.(df)/3),
                           vmin=0)

    # mutational sensitivity map
    @assert size(df, 2) == 2 # make sure data is formatted properly
    fig, ax = subplots(figsize=(5,4))
    formatFig!(fig)
    fig.subplots_adjust(left = 0.1,bottom=0.1, top=0.9, right=0.90)
    fs = 16
    orientNetwork!(Q)
    df = vec(mean( abs.(df), dims=2))
    df_matrix = bondDataVector2Matrix(df, Q.A)
    plotColoredBonds!(ax, Q.structs[1], Q.A, df_matrix, vmax, cmap="Blues")

    ax.axis("equal")
    ax.axis("off")
    makeColorbar!(fig, cm=cmap, low=0.0, high=vmax, x=.92, y=.20, wid=.04, len=.6)
    
    fig.text(.99, .5, "Δ $assay", rotation=270, fontsize=16,
             horizontalalignment="center", verticalalignment="center")
    return fig, "bondSensMap_" * assay
end





function genColoredStructFig(Q::Network,
                             df::Vector{Float64},
                             assay::String;
                             cmap="Blues",
                             vmax=maximum(abs.(df)),
                             vmin=minimum(abs.(df)))

    # mutational sensitivity map
    fig, ax = subplots(figsize=(5,4))
    formatFig!(fig)
    fig.subplots_adjust(left = 0.1, bottom=0.1, top=0.9, right=0.80)
    fs = 16
    orientNetwork!(Q)
    plotColoredStruct!(ax, Q.structs[1], Q.A, df, cmap=cmap, vmax=vmax)

    ax.axis("equal")
    ax.axis("off")
    makeColorbar!(fig, cm=cmap, low=vmin, high=vmax, x=.85, y=.20, wid=.04, len=.6)
    
    fig.text(.95, .5, assay, rotation=270, fontsize=16,
             horizontalalignment="center", verticalalignment="center")
    return fig
end





function genDMSFig(df::Matrix{Float64},
                   assay::String;
                   colorTop=maximum(abs.(df)))
    # plot deep mutational scan arrays

    numNodes, numTypes = size(df)
    figx = numNodes*.25 + 1
    figy = numTypes*.25 + 4

    fig, ax = subplots(figsize=(figx,0.6figy))
    formatFig!(fig)
    fig.subplots_adjust(left=0.1, bottom=0.1, top=0.9, right=0.85, hspace=0.1)
    top = .05; fs = 16    
    h = [.45, -2.1]
    plotDMS!(ax, df, colorTop)
    ax.set_ylabel("Node Type")
    ax.set_xlabel("Position")

    ax.text(1.12, h[1], "Δ $assay", transform=ax.transAxes,
           rotation=270, fontsize=14, horizontalalignment="center",
          verticalalignment="center")
    makeColorbar!(fig, cm="bwr", low=-colorTop, high=colorTop, x=.855, y=.3, wid=.02, len=.39)
    return fig, "DMS_"*assay
end


function genConfChangeDMSFig(df::Matrix,
                    assay::String;
                    colorTop=maximum(abs.(df)))
    # plot confromational change deep mutational scan
    numNodes, numTypes = size(df)
    figx = numNodes*.25 + 1
    figy = numTypes*.25 + 4

    fig, ax = subplots(figsize=(figx,0.6figy))
    formatFig!(fig)
    fig.subplots_adjust(left=0.1, bottom=0.1, top=0.9, right=0.85, hspace=0.1)
    top = .05; fs = 16    
    h = [.45, -2.1]
    plotDMS!(ax, df, colorTop)
    ax.set_ylabel("Node Type")
    ax.set_xlabel("Position")
    ax.text(1.08, h[1], "Conf. Change", transform=ax.transAxes,
           rotation=270, fontsize=14, horizontalalignment="center",
          verticalalignment="center")
    makeColorbar!(fig, cm="bwr", low=-colorTop, high=colorTop, x=.855, y=.3, wid=.02, len=.39)
    return fig, "CC_DMS_"*assay
end; export genConfChangeDMSFig



function genMutSensDistFig(df::Matrix{Float64},
                           assay::String)
    
    fig, ax = subplots(figsize=(6,3))
    formatFig!(fig)
    fig.subplots_adjust(left=0.1, bottom=0.2, top=0.9, right=0.9, hspace=.15)

    # x lims 
    m, M = extrema(df)
    xrange = [m,M]
    ymax = Float64[]
    bins = linspace(m,M,25)
    plotMutSensHist!(ax, df, bins)
    ax.axvline(x=0,color="r", zorder=1, label="_nolegend_")
    push!(ymax,ax.get_ylim()[2])
    ax.set_xlim(xrange)
    for i in 1:2; ax.set_ylim([0,maximum(ymax)*1.1]);end
    ax.set_xlabel("Δ $assay") 

   # Mixture = pyimport("sklearn.mixture")
   # data = reshape(filter(x->x!=0, df[:]), :,1)
   # gmm = Mixture.GaussianMixture(n_components=2,covariance_type="full").fit(data)
   # gaussian(x,μ,σ) = (1/(σ*sqrt(2*pi))) * exp(-0.5*((x-μ)/σ)^2)
   # μ1, μ2 = gmm.means_
   # σ1, σ2 = sqrt.(gmm.covariances_)
   # w1, w2 = gmm.weights_
   # dgparams = μ1, μ2, σ1, σ2, w1, w2
   # x = LinRange(m,M,1000)
   # y = doubleGaussian.(x, μ1, μ2, σ1, σ2, w1, w2)
   # y1 = w1*gaussian.(x,μ1,σ1)
   # y2 = w2*gaussian.(x,μ2,σ2)
   # ax.plot(x,y1)
   # ax.plot(x,y2)
   # ax.plot(x,y)
    return fig, "mutSensDist_"*assay
end


function genMutSensRDFFig(df::Matrix{Float64},
                          adjMat::Matrix{Float64},
                          activeSiteNodes::Vector{<:Number},
                          assay::String)

    LSM = computeLeastStepMatrix(adjMat)
    steps = computeStepsFromAct(LSM, actNodes=activeSiteNodes)
    ps = 120
    al = .2

    fig, ax = subplots(figsize=(5,5))
    formatFig!(fig)
    fig.subplots_adjust(left=0.2, bottom=0.1, top=0.9, right=0.9)
    fig.subplots_adjust(wspace=.1)
    x = repeat(steps,1,size(df,2))
    y = abs.(df)
    ax.scatter(x,y,ps,color="k", alpha=al, edgecolor="")
    ax.set_yscale("log")
    ax.set_ylim([8e-6, 5e-1])
    ax.set_axisbelow(true)
    ax.grid(linestyle="--")
    ax.set_xlabel("Distance from active site")
    ax.set_ylabel("Δ $assay")
    return fig, "mutSensRDF_" * assay
end



function genStabFitDMSFig(stab::Vector{Float64}, func::Vector{Float64},
                          stabWT::Float64, funcWT::Float64;
                          xlabel::AbstractString="Stability",
                          ylabel::AbstractString="Binding")
    fig, ax = subplots(figsize=(5,3.5))
    fig.subplots_adjust(left = 0.25, bottom=0.2, top=0.9, right=0.9)
    formatFig!(fig)
    scatter(stab, func, 20, color="k")
    scatter([stabWT], [funcWT], 20, color="r")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_axisbelow(true)
    ax.grid(linestyle="--")
    return fig
end


function genStabBindFig(stab::Vector{Float64}, bind::Vector{Float64},
                          stabWT::Float64, bindWT::Float64;
                          unfoldedEnergy::Float64=0.5)
    # this function takes fitness values as arguments not biophysical values
    return genStabFitDMSFig(-stab .- unfoldedEnergy,
                            -bind,
                            -stabWT - unfoldedEnergy,
                            -bindWT,
                            xlabel=L"ΔE_{fold}",
                            ylabel=L"ΔE_{bind}")
end


function genStabBindScatterHistFig(stab::Vector{Float64}, bind::Vector{Float64},
                                   stabWT::Float64, bindWT::Float64;
                                   unfoldedEnergy::Float64=0.5)
    fig = scatterHist(-stab .- unfoldedEnergy, -bind, nbins=30)
    formatFig!(fig)
    ax = fig.axes[1]
    ax.scatter([-stabWT - unfoldedEnergy], [-bindWT], 20, c="r")
    ax.axvline(x=(-stabWT - unfoldedEnergy), c="r", lw=.5, zorder=0)
    ax.axhline(y=-bindWT, c="r", lw=.5, zorder=0)
    ax.set_xlabel(L"ΔE_{fold}")
    ax.set_ylabel(L"ΔE_{bind}")
    ax.text(.75,.94, L"E_{unfolded} = "*"$unfoldedEnergy", transform=ax.transAxes)
    return fig, "stab_bind"
end


function genEvolFig(evoTraces::Matrix{Float64}, assay::String)
    # fitness traces in the columns of evoTraces
    fig, ax = subplots(figsize=(8,5))
    formatFig!(fig)
    fig.subplots_adjust(left = 0.2,bottom=0.1, top=0.9, right=0.9)

    famSize = size(evoTraces, 2)
    alpha = sqrt(famSize) / famSize
    ax.plot(evoTraces, c="k", alpha = alpha, label="_nolegend_") 
   

    avgline = vec(mean(evoTraces, dims=2))
    ax.plot(avgline, linewidth=4, color="k")
    ax.set_xlabel("Evolutionary time")
    ax.set_ylabel(assay)
    ax.legend(["average"])
    return fig
end



function genEnergyVsConfCoorFig(Q::Network)
    SM, RLM, EM = Q.pheno
    orientNetwork!(Q)
    dxy = Q.structs[5] .- Q.structs[1]
    xy = similar(dxy)
    N = 100
    #x = linspace(-.9, 1.9, N)
    x = linspace(-.4, 1.4, N)
    
    # solvent - solvent
    addLig!(Q, 1, 1) # put solvent at active site
    addLig!(Q, 1, 2) # put solvent at allosteric site
    energies_0 = zeros(N)
    for (i,a) in enumerate(x)
        xy .= Q.structs[1] .+ a .* dxy 
        energies_0[i] = computeEnergy(xy, SM, RLM, EM)
    end

    # ligand ligand
    addLig!(Q, 2, 1) # put ligand at active site
    addLig!(Q, 2, 2) # put ligand at allosteric site
    energies_1 = zeros(N)
    for (i,a) in enumerate(x)
        xy .= Q.structs[1] .+ a .* dxy 
        energies_1[i] = computeEnergy(xy, SM, RLM, EM)
    end

    fig, ax = subplots(figsize=(2.6,2))
    formatFig!(fig)
    ax.plot(x, energies_0, c=[0,0,1], linewidth=3)
    ax.plot(x, energies_1, c=[.1,.65,.1], linewidth=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([])
    #ax.set_xlabel("Conf. Coordinate
    #ax.set_ylabel("energy")
    return fig, "energy_vs_confCoor"
end


function genEnergyVsConfCoorFig2(Q::Network)
    SM, RLM, EM = Q.pheno
    orientNetwork!(Q)
    dxy = Q.structs[5] .- Q.structs[1]
    xy = similar(dxy)
    N = 100
    #x = linspace(-.9, 1.9, N)
    x = LinRange(-.4, 1.4, N)
    
    # solvent - solvent
    removeAllLig!(Q)
    energies_0 = zeros(N)
    for (i,a) in enumerate(x)
        xy .= Q.structs[1] .+ a .* dxy 
        energies_0[i] = computeEnergy(xy, SM, RLM, EM)
    end


    fig, ax = subplots(figsize=(2.6,2))
    formatFig!(fig)
    ax.plot(x, energies_0, c=[0,0,1], linewidth=3)
    #ax.plot(x, energies_1, c=[.1,.65,.1], linewidth=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([])
    ax.set_xlabel("Conf. Coordinate")
    ax.set_ylabel("energy")
    return fig, "energy_vs_confCoor"
end




function genEnergySpectrumFig(energies::Vector{<:Real})
    N = 40
    if length(energies) > N
        energies = energies[1:N]
    end
    
    sort!(energies)
    fig, ax = subplots()
    formatFig!(fig)
    ax.scatter(1:length(energies), energies, c="k", marker="_")
    ax.set_ylabel("energy")
    ax.set_xlabel("state index")
    ax.set_axisbelow(true)
    ax.grid(linestyle="--")
    return fig
end




function genLigEnergySliceFig(ligLengths::AbstractVector{<:Number},
                              energies_solAtAllo::AbstractVector{<:Number},
                              energies_ligAtAllo::AbstractVector{<:Number},
                              Q::Network)

    fig, ax = subplots(figsize=(4,3))
    formatFig!(fig)
    fig.subplots_adjust(left=0.2, bottom=0.2, top=0.9, right=0.9)


    ax.plot(ligLengths, energies_solAtAllo, c="k", linewidth=2 )
    ax.plot(ligLengths, energies_ligAtAllo, "--", c="k", linewidth=2 )
    ax.set_xlabel("ℓ", fontsize=16)
    ax.set_ylabel("energy", fontsize=16)
    ax.set_axisbelow(true)
    ax.grid(linestyle="--")
    ax.axvline(x=Q.ligs[1].bonds[1].l, color=[.5,.5,1.],ls="--",lw=2, zorder=1, label="_nolegend_")
    ax.axvline(x=Q.ligs[2].bonds[1].l, color=[.5,.8,.5],ls="--",lw=2, zorder=1, label="_nolegend_")
    ax.legend(["Solvent","Ligand"])
    return fig, "LigSlice"
end

function genLigScapeFig(energies::Matrix{Float64},
                            ligLengths::AbstractVector{<:Number},
                            ls::Number, lr::Number )
    fig, ax = subplots()
    fig.subplots_adjust(left=0.1, bottom=0.2, top=0.9, right=0.9)
    cp = ax.contourf(ligLengths, ligLengths, energies, cmap="YlGn_r", antialiased=true)
    clb = fig.colorbar(cp)
    #clb.set_label("Energy", fontsize=24, rotation=-90,labelpad=30)
    #ax.set_xlabel(L"\ell_{allosteric}", fontsize=26)
    #ax.set_ylabel(L"\ell_{active}", fontsize=26)
    ligs = [ls,lr]
    x = repeat(ligs,2)
    y = repeat(ligs',2)[:]
    ax.scatter(x,y, 100, c="k", marker="+")
   # ax.axis("equal")
    ax.set_xlim([minimum(ligLengths), maximum(ligLengths)])
    ax.set_ylim([minimum(ligLengths), maximum(ligLengths)])
    ax.set_aspect("equal")
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=14)
    #ax.set_xticklabels(fontdict=Dict("fontsize"=>14))
    return fig, "ligScape"
end

function genLigScapeFig(energies::Matrix{Float64},
                            ligLengths::AbstractVector{<:Number},
                            Q::Network)
    ls = Q.ligs[1].bonds[1].l
    lr = Q.ligs[2].bonds[1].l
    return genLigScapeFig(energies, ligLengths, ls, lr)
end



function genLigScapeFig(energies::Matrix{Float64},
                            structs::Matrix{Matrix{Float64}},
                            ligLengths::AbstractVector{<:Number},
                            Q::Network)
    # this is the version that 
    stateLables = labelLigScape(Q, structs)
    xs, ys = outlineLigScape(stateLables, ligLengths);
    ls = Q.ligs[1].bonds[1].l
    lr = Q.ligs[2].bonds[1].l
    fig, name = genLigScapeFig(energies, ligLengths, ls, lr)
    ax = fig.axes[1]
    ax.plot.(xs, ys, "k:")  # plot outlines
    return fig, "ligScape"
end


#function genLigScapeAlloLabelStatesFig(labels::Matrix{Int},
#                                       ligLengths::AbstractVector{<:Number})
#    fig, ax = subplots(figsize=(5,5))
#    l = ligLengths[1]
#    L = ligLengths[end]
#    ax.imshow(labels, cmap="nipy_spectral", origin="lower",
#              extent=(l,L,l,L))
#    ax.set_xlabel("Allosteric site ligand length", fontsize=16)
#    ax.set_ylabel("Active site ligand length", fontsize=16)
#    ls = 1.3 
#    lr = 1.9 
#    ligs = [ls,lr]
#    x = repeat(ligs,2)
#    y = repeat(ligs',2)[:]
#    ax.scatter(x,y, 100, c="k", marker="+")
#    ax.axis("equal")
#   #ax.set_xlim([minimum(ligLengths), maximum(ligLengths)])
#   # ax.set_ylim([minimum(ligLengths), maximum(ligLengths)])
#    return fig
#end

function genEpiMatrixFig(epiMat::Matrix{<:Number}, assay::String)
    N = size(epiMat,1)
    fig, ax = subplots()
    x = 0:5:N
    ax.set_xticks(x)
    ax.set_xticklabels(x)
    ax.set_yticks(x)
    ax.set_yticklabels(reverse(x) .+ 5)
    im = ax.imshow(epiMat, extent=[1,N,1,N])
    clb = fig.colorbar(im,fraction=0.046, pad=0.03)
    mn, mx = extrema(epiMat)
    clb.set_ticks([mn, mx])
    clb.set_ticklabels([@sprintf("%3.2f", mn),@sprintf("%3.2f",mx)])
    clb.set_label("ΔΔ $assay", fontsize=14, rotation=-90)
    return fig
end

function genEpiHistFig(epi::Vector{Matrix{Float64}},
                       assay::String)
    # epi is the array structure that stores all the double mutation
    # epistatic terms.
    
    # format data
    N = length(epi); M = length(epi[1])
    E = zeros(N*M)
    for i in 1:N, j in 1:M
        E[M*(i-1)+j] = epi[i][j]
    end

    # plot histogram
    fig, ax = subplots(figsize=(6,4))
    formatFig!(fig)
    fig.subplots_adjust(left=0.1, bottom=0.2, top=0.9, right=0.9)
    mn, mx = extrema(E)
    bns = linspace(mn,mx,100)
    ax.hist(E,color="k", rwidth=1, bins=bns, zorder=2)
    ax.set_axisbelow(true)
    ax.grid(linestyle="--")
    ax.set_ylabel("Counts")
    ax.set_title("Epistatic terms",fontsize=16)
    ax.set_xlabel("ΔΔ $assay")
    ax.set_yscale("log")
    return fig
end

function genEpiMatSpectrumFig(epiMat::Matrix{Float64}, assay::String)

    evals = eigvals(epiMat)
        
    # plot histogram
    fig, ax = subplots(figsize=(6,4))
    formatFig!(fig)
    fig.subplots_adjust(left=0.1, bottom=0.2, top=0.9, right=0.9)
    mn, mx = extrema(evals)
    bns = linspace(mn,mx,50)
   
    ax.hist(evals, color="k", rwidth=1, bins=bns, zorder=2)
    ax.set_axisbelow(true)
    ax.grid(linestyle="--")
    ax.set_ylabel("Counts")
    ax.set_title("Epistatic matrix eigen spectrum, $assay",fontsize=16)
    ax.set_xlabel("λ", fontsize=16)
    return fig
end


############################
# SCA figures ##############
############################

    
function genSeqSimFig(M::Matrix{<:Number})
    fig, ax = subplots()
    fig.subplots_adjust(left=0.1, right=0.8, bottom=0.1, top=0.9)
    im = ax.imshow(M, vmin=0, vmax=1)
    clb = fig.colorbar(im,fraction=0.046, pad=0.03)
    clb.set_label("Hamming Distance", fontsize=14, rotation=-90, labelpad=20)
    ax.set_title("Sequence Similarity", fontsize = 18)
    return fig
end

function genFitnessesFig(fits::Matrix{<:Number},
                         assay=nothing)
    numEvoSteps, famSize = size(fits)
    fig, ax = subplots()
    alpha = 1/ sqrt(famSize)
    ax.plot(fits, c="k", alpha=alpha)
    ax.set_xlabel("Evolutionary time", fontsize=16)
    if assay == nothing
        ax.set_ylabel("Fitness", fontsize=16)
    elseif isa(assay, String)
        ax.set_ylabel(assay, fontsize=16)
    else
        error("Wrong assay variable type")
    end
    return fig
end

function genConsMapFig(Q::Network, cons::Vector)
    fig = genColoredStructFig(Q, cons, "Conservation", vmax=1.6, cmap="Purples")
   # fig.subplots_adjust(left=0.1, right=0.7, bottom=0.1, top=0.9)
    return fig
end

function genPosMatrixFig(M::Matrix{<:Number}, assay::String; vmax=maximum(M))
    N = size(M,1)
    fig, ax = subplots()
    fig.subplots_adjust(left=0.1, right=0.8, bottom=0.1, top=0.9)    
    x = 0:5:N
    ax.set_xticks(x)
    ax.set_xticklabels(x)
    ax.set_yticks(x)
    ax.set_yticklabels(reverse(x) .+ 5)
    im = ax.imshow(M, extent=[1,N,1,N], vmax=vmax)
    clb = fig.colorbar(im,fraction=0.046, pad=0.03, format="%.0e")
    mx = vmax
    mn = minimum(M)
    xs = linspace(mn, mx, 10)
    clb.set_label("$assay",labelpad=15, fontsize=14, rotation=-90)
    return fig
end

#######################################


function genStrainVsSoftModeFig(Q)
    
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
    overlaps = (E.vectors' * dx)[4:end]
    values = E.values[4:end]
    val, rank = findmax(abs.(overlaps))
    rank = 4

    m = -E.vectors[:,rank]
    dxy = Array(reshape(dx, 2, :)')
    mm = Array(reshape(m, 2, :)')
    dot(dx, m) < 0 && (mm .*= -1)


    # map figure
    fig, ax = subplots(figsize=(6,6))
    fig.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)
    plotNetwork!(ax, xy1, Q.A, bondColor=[.7,.7,.7])
    ax.axis("off")
    ax.axis("equal")
    ax.quiver(xy1[:,1], xy1[:,2], dxy[:,1], dxy[:,2], zorder=2, scale=2 )

    ax.quiver(xy1[:,1], xy1[:,2], mm[:,1], mm[:,2], color="r", zorder=2.2, scale=2.5)

    xlim = ax.get_xlim()
    ylim= ax.get_ylim()
    xlim = xlim .+ [-0.4, 0.4]
    ylim = ylim .+ [-0.4, 0.4]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


    # plot figure
    N = length(overlaps)
    fig2, ax2 = subplots(figsize=(4,3))
    fig2.subplots_adjust(left=0.22, right=0.9, bottom=0.22, top=0.9)
    formatFig!(fig)
    #ax2.plot(values, zeros(N), c="k")
    ax2.vlines(x=values, ymin=0, ymax=abs.(overlaps), color="k")
    ax2.plot(values, abs.(overlaps), ".k", ms=10)
    #ax2.stem(1:length(overlaps),-overlaps,  markerfmt=".k", linefmt="k", basefmt="k")
    ax2.tick_params(axis="both", which="major", labelsize=16)
    ax2.set_yticks([0,0.5,1])
    ax2.set_yticklabels([0,.5,1])
    ax2.set_ylim([0,1.05])
    ax2.set_xlim([5e-3,2e1])
    ax2.set_xscale("log")
    ax2.set_xticks([1e-2,1e-1,1e0, 1e1])
    ax2.set_xlabel(L"\lambda_k", fontsize=16)
    ax2.set_ylabel(L"q_k", fontsize=16)
    return fig, "softModeAlloModeMap", fig2, "modeSpectrum"
end



function genRegistryFig(reg::Vector{Matrix{Float64}}, A::Matrix;
                        x=3,
                        y=5)
    fig, ax = subplots(x,y, figsize = (y*5, x*4))
    num = x * y
    for i in 1:num
        plotNetwork!(ax[i], reg[i], A, bondColor="grey")
        ax[i].axis("off")
        ax[i].axis("equal")
    end
    return fig
end














