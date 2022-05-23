
#############################################
#######  PLOTTING FUNCTIONS  ################
#############################################

function plotNetwork(xy, A)
    fig, ax = subplots()
    plotNetwork!(ax, xy, A)
    return fig
end

function plotNetwork(Q::Network)
    return plotNetwork(Q.structs[1], Q.A)
end


function plotNetwork!(ax, xy::Matrix, adjacencyMatrix::Matrix; bondColor=clr("blue"),
                     bondWidth=3, nodeSize=30)
    # Plots network, modifies figure axis object in place
    # uses PyPlot, LinearAlgebra
    x = xy[:,1]
    y = xy[:,2]
    N = size(adjacencyMatrix,1)
    A = LowerTriangular(adjacencyMatrix)
    
    for i =1:N
        whichj = findall(A[i,:].>0)
        for j in whichj
            
            if A[i,j] == 1
                ax.plot([x[i],x[j]], [y[i],y[j]], lw=bondWidth, c=bondColor, zorder=1)
            elseif A[i,j] < 1
                ax.plot([x[i],x[j]], [y[i],y[j]], ":", lw=bondWidth, c=bondColor, zorder=1)
            end
        end
    end
    ax.scatter(x,y, nodeSize, c="black",zorder=2)
    return nothing
end


function plotComplex(Q::Network, xy::Matrix;
                     actSite::CartesianIndex{2}=CartesianIndex(3, 5),
                     alloSite::CartesianIndex{2}=CartesianIndex(1, 7),
                     bondWidth=3, nodeSize=30,
                     xlim=[minimum(xy[:,1])-1, maximum(xy[:,1])+1], 
                     ylim=[minimum(xy[:,2])-1, maximum(xy[:,2])+1])

    A = Q.A
    numNodes = size(Q.A, 1)
    numLigs = length(Q.ligs)
    actSite = Q.sites[1]
    alloSite = Q.sites[2]
    colors=[[0,0,1], [0.1, .65 ,0.1]]

    fig, ax = subplots(figsize=(5,4))
    formatFig!(fig); fs = 16
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)
    plotNetwork!(ax, xy, A, bondColor=[.5,.5,.5], bondWidth=bondWidth, nodeSize=nodeSize)
    ax.axis("off")
    ax.axis("equal")

    # plot active site ligand
    for (j,loc) in enumerate(actSite)
        ax.plot( [xy[loc[1],1], xy[loc[2],1]],
                   [xy[loc[1],2], xy[loc[2],2]],
                   color=colors[1], linewidth=bondWidth+2, zorder=1 )
    end

    # plot active site ligand
    for (j,loc) in enumerate(alloSite)
        ax.plot( [xy[loc[1],1], xy[loc[2],1]],
                   [xy[loc[1],2], xy[loc[2],2]],
                   color=colors[1], linewidth=bondWidth+2, zorder=1 )
    end
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return fig
end



function plotComplex!(ax, xy, adjacencyMatrix, numNetNodes;
                     netBondColor=clr("blue"), ligBondColor=clr("green"),
                     ligNetBondColor=clr("red"))
    # Plots Ligand Network complex, returns figure object
    # uses PyPlot, LinearAlgebra

    x = xy[:,1]
    y = xy[:,2]
    N = size(adjacencyMatrix,1)
    A = LowerTriangular(adjacencyMatrix)

    for i =1:N
        whichj = findall(A[i,:].>0)
        for j in whichj
            if i <= numNetNodes && j <= numNetNodes
                ax.plot([x[i],x[j]], [y[i],y[j]], lw=3, c=netBondColor, zorder=1)
            elseif i > numNetNodes && j <= numNetNodes
                ax.plot([x[i],x[j]], [y[i],y[j]], lw=2, c=ligNetBondColor, zorder=1)
            elseif i > numNetNodes &&  j > numNetNodes
                ax.plot([x[i],x[j]], [y[i],y[j]], lw=3, c=ligBondColor, zorder=1)
            end
        end
    end
    xMax = maximum(x)
    yMax = maximum(y)
    xMin = minimum(x)
    yMin = minimum(x)
    ax.set_xlim(round(xMin-.5),round(xMax+.5))
    ax.set_ylim(round(yMin-.5),round(yMax+.5))
    ax.scatter(x,y,65,color=clr("black"),zorder=2)
    ax.axis("off")
    return nothing
end


function plotColoredStruct!(ax,
                            xy::Matrix{Float64},
                            adjacencyMatrix::Matrix{Float64},
                            score::Vector{<:Number};
                            cmap::String="Blues",
                            vmax=0.05,
                            vmin=0)
    
    # Plots network, returns figure object
    x = xy[:,1]
    y = xy[:,2]
    numNodes = size(adjacencyMatrix,1)
    A = LowerTriangular(adjacencyMatrix)
    if vmax < minimum(score)
        vmax = minimum(score)
    end
    
    # plot Bonds
    for i =1:numNodes
        whichj = findall(A[i,:].>0)
        for j in whichj
            ax.plot([x[i],x[j]], [y[i],y[j]], lw=5, c=clr("grey"), zorder=1)
        end
    end
    ax.scatter(x,y,s=200, c=clr("bblack"), zorder=3)
    ax.scatter(x,y,s=150, c=score, cmap=cmap, vmax=vmax, vmin=vmin, zorder=4)
    return nothing
end



function plotDMS!(ax, df, top)
    numNodes = size(df,1)
    numTypes = 5 #size(df,2)
    df_norm = zeros(size(df))
    for i in eachindex(df)
        if df[i] >= top
            df_norm[i] = 100 
        elseif df[i] < top && df[i] > -top
            df_norm[i] = Int(round((50*df[i]/top) + 50))
        elseif df[i] <= -top
            df_norm[i] = 0
        end
    end
    
    ax.imshow(df', cmap="bwr", vmin=-top, vmax=top)
    ax.tick_params(length=0)
    ax.set_xticks(0:numNodes-1)
    ax.set_xticklabels(1:numNodes, fontsize=9)
    ax.set_ylim(-.5, numTypes-.5)
    ax.set_yticks(0:numTypes-1)
    ax.set_yticklabels(numTypes:-1:1)
    return nothing
end




function plotMutSensHist!(ax, df, bns)

    ax.axvline(x=0,color="r", zorder=1, label="_nolegend_")
    #ax.hist(filter(x->x!=0,df[:]), color="k", rwidth=.94, bins=bns, zorder=2)
    ax.hist(filter(x->x!=0,df[:]), color="k", rwidth=.94, bins=bns, zorder=2, density=true)
    ax.set_axisbelow(true)
    ax.grid(linestyle="--")
    ax.set_ylabel("Counts")
    return nothing
end


function plotColoredBonds!(ax,
                           xy::Matrix{Float64},
                           A::Matrix{Float64},
                           bondDataMatrix::Matrix{Float64},
                           top::Number;
                           cmap::String="Greens",
                           numShades::Int=100)
    # plots bond strains from ligand binding on network 
    x = xy[:,1]
    y = xy[:,2]
    N = size(A,1)
    A = LowerTriangular(A)
    bondDataMatrix = abs.(bondDataMatrix)
    cmp = get_cmap(cmap,numShades) # a pyplot function

    for i =1:N
        whichj = findall(A[i,:].>0)
        for j in whichj
            if bondDataMatrix[i,j] <= top
                colorIndex = Int(round((numShades*bondDataMatrix[i,j]/top))+1)
            elseif bondDataMatrix[i,j] > top
                colorIndex = numShades
            end
            ax.plot([x[i],x[j]], [y[i],y[j]], lw=6, c=clr("grey"), zorder=1)
            ax.plot([x[i],x[j]], [y[i],y[j]], lw=5, c=cmp(colorIndex), zorder=1)
        end
    end
    ax.scatter(x,y,s=100,c=clr("ggrey"),zorder=2)
    return top
end



function plotDecay!(ax, x, y; color = "r")
    ax.scatter(x, y, s=125, c=color,alpha=.3, edgecolor="none")
    ax.set_xscale("log")
    ax.set_yscale("log")
end


function plotEnergyTable!(ax, E_table; withLines=true)
    labels = ["solvent", "right", "wrong"]
    x = [1,2,3]
    n = size(E_table,1)
    ms = 1000

    ax.scatter(repeat([1.],n), E_table[:,1], color="b", marker="_", ms, zorder=2)
    ax.scatter(repeat([2.],n), E_table[:,2], color="g", marker="_", ms, zorder=2)
    ax.scatter(repeat([3.],n), E_table[:,3], color="r", marker="_", ms, zorder=2)

    if withLines
        for i in 1:n
            ax.plot([1.21,1.79], E_table[i,1:2],linewidth=1, c=[.6,.6,.6], zorder=1)
            ax.plot([2.21,2.79], E_table[i,2:3],linewidth=1, c=[.6,.6,.6], zorder=1)
        end
    end
    ax.set_xlim([.7,3.3])
    ax.set_xticks(x)
    ax.set_ylabel("Energy",fontsize=16)
    ax.set_xticklabels(labels, fontsize=16)
end



function plotConfChangeFreqs!(ax, counts0, counts1)
    # fractions = [ fraction two state before evolution
    #             , fraction two state after evolution]
    y = [counts0, counts1]
    ax.bar([1,2], y, color = "k")
    ax.set_xticks([1,2])
    ax.set_xticklabels(["Before","After"], fontsize = 16)
    ax.hlines(y=200, xmin=.6, xmax=1.4)
    #ax.hlines(y=200, xmin=.6, xmax=1.4)
end



function plotBarChart(data::Array{Array{Float64,1},1},
                      sampleSizes::Array{<:Number,1},
                      titles::Array{String,1},
                      ylbl::String)

    numEns = length(data)
    figx = 3*numEns
    figy = 3

    fig, ax = subplots(1, numEns, figsize=(figx, figy))
    numEns == 1 && (ax = reshape([ax],1,1))
    formatFig!(fig)
    fig.subplots_adjust(left = 0.15, bottom=0.1, top=0.9, right=0.9, wspace=0.3, hspace=0.05)

    for i in 1:numEns
        y = data[i]
        ax[i].bar([1,2], y, color = "k")
        ax[i].set_xticks([1,2])
        ax[i].set_xticklabels(["Before","After"], fontsize = 16)
        ax[i].set_ylim([0, 1.0])
        ax[i].text(0.08,0.92,"N = $(Int(sampleSizes[i]))", transform=ax[i].transAxes)
        ax[i].set_title("Selected for "*titles[i], fontsize=14)
    end
    ax[1].set_ylabel(ylbl)


    setLimSame!(ax,"y")
    return fig, ax
end


function plotScatter!(ax,
                      x::Matrix{<:Number},
                      y::Matrix{<:Number})

    ps = 45; a = .35
    ax.scatter(x[:,1], y[:,1], ps, color="grey", alpha=a, edgecolors="")
    ax.scatter(x[:,2], y[:,2], ps, color="g", alpha=a, edgecolors="")
    ax.legend(["Before","After"], loc=2)
    ax.set_axisbelow(true)
    ax.grid(linestyle="--")
    return nothing
end


function plotCorrelations(datax::Vector{<:Matrix{<:Real}},
                          datay::Vector{<:Matrix{<:Real}},
                          titles::Array{String,1},
                          titleFontSize=16)

    @assert length(titles) == length(datax)
    numEns = length(datax)
    figx = 5*numEns; figy = 4
    fig, ax = subplots(1, numEns, figsize=(figx, figy))
    numEns==1 && (ax=reshape([ax],1,1))
    formatFig!(fig)
    fig.subplots_adjust(left = 0.05, bottom=0.2, top=0.9, right=0.95, wspace=0.15, hspace=0.1)
    for i in 1:numEns
        plotScatter!(ax[i], datax[i], datay[i])
        ax[i].set_title("Selected for "*titles[i], fontsize=titleFontSize)
    end
    return fig, ax
end



function plotCorrelations(datax::Vector{<:Matrix{<:Real}},
                          datay::Matrix{<:Matrix{<:Real}},
                          titles::Array{String,1};
                          titleFontSize=16)
    @assert length(titles) == length(datax)
    numAssay, numEns = size(datay)
    figx = 5*numEns
    figy = 3*numAssay
    fig, ax = subplots(numAssay, numEns, figsize=(figx, figy))
    (numAssay==1 && numEns==1) && (ax=reshape([ax],1,1))
    formatFig!(fig)
    fig.subplots_adjust(left = 0.15, bottom=0.1, top=0.9, right=0.9, wspace=0.15, hspace=0.1)

    for i in 1:numAssay
        for j in 1:numEns
            plotScatter!(ax[i,j], datax[j], datay[i,j])
            if i==1; ax[i,j].set_title("Selected for "*titles[j], fontsize=titleFontSize); end
            if i<numAssay; ax[i,j].set_xticklabels([]); end
            if j>1; ax[i,j].set_yticklabels([]); end
        end
    end
    return fig, ax
end




function scatterHist(x, y; nbins=30, equalAxes=false)

    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    hist_height = 0.15
    spacing = 0.0
    buffer = 0.1
    barwid = 0.9

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, hist_height]
    rect_histy = [left + width + spacing, bottom, hist_height, height]

    dx = 6; dy = 5
    fig = figure(figsize=(dx, dy))
    ax = fig.add_axes(rect_scatter)
    ax.scatter(x, y, 10, c="k")
    ax.set_axisbelow(true)
    ax.grid(linestyle="--")

    xl = extrema(x)
    yl = extrema(y)
    xl = (xl[1] - buffer * (xl[2]-xl[1]), xl[2] + buffer * (xl[2]-xl[1]))
    yl = (yl[1] - buffer * (yl[2]-yl[1]), yl[2] + buffer * (yl[2]-yl[1]))
    Dx = xl[2] - xl[1]
    Dy = yl[2] - yl[1]
    if (Dy / dy) > (Dx / dx)
        α = 0.5 * ( (dx/dy)*(Dy) - Dx )
        new_xl = (xl[1] - α, xl[2] + α) 
        new_yl = yl
    else
        α = 0.5 * ( (dy/dx)*(Dx) - Dy )
        new_xl = xl 
        new_yl = (yl[1] - α, yl[2] + α)
    end
    ax_histx = fig.add_axes(rect_histx)
    ax_histy = fig.add_axes(rect_histy)
    ax_histx.tick_params(axis="x", labelbottom=false)
    ax_histx.set_yticks([])
    ax_histx.set_yticklabels([])
    ax_histy.tick_params(axis="y", labelleft=false)
    ax_histy.set_xticks([])
    ax_histy.set_xticklabels([])
    ax_histx.spines["right"].set_visible(false)
    ax_histx.spines["top"].set_visible(false)
    ax_histx.spines["left"].set_visible(false)
    ax_histy.spines["bottom"].set_visible(false)
    ax_histy.spines["top"].set_visible(false)
    ax_histy.spines["right"].set_visible(false)
    xbins = linspace(new_xl[1],new_xl[2],nbins)
    ybins = linspace(new_yl[1], new_yl[2],nbins)
    ax_histx.hist(x, bins=xbins,color="k", rwidth=barwid)
    ax_histy.hist(y, bins=ybins,color="k", rwidth=barwid, orientation="horizontal")
    ax.set_xlim(new_xl)
    ax.set_ylim(new_yl)
    ax_histx.set_xlim(new_xl)
    ax_histy.set_ylim(new_yl)
    return fig
end


function computeWilsonCI(n::Int, p::Number, α::Number)
    # compute the Wilson Confidence interval for  Binomial data 
    z = cquantile(Normal(), α/2)
    d = (1 + z^2/n)
    center = (p + z^2/(2n)) / d
    spread = sqrt( p*(1-p)/n + z^2/(4n^2) )  / d
    lower = center - spread
    upper = center + spread
    return lower, upper
end

function wilsonCIHelper(df::Vector, α::Number, ind::Int)

    # this function takes the raw phenotypes data and 
    # returns the confidence intervals curves.
    CI = zeros(length(df), 2)
    for i in 1:length(df)
        data = df[i][:,ind]
        CI[i,:] .= computeWilsonCI(length(data), mean(data), α)
    end
    return CI
end

function computeWaldCI(n::Int, μ::Number, σ::Number, α::Number)
    # n = size of sample
    # μ = mean of sample
    # σ = std of sample
    # α = confidence
    z = cquantile(Normal(), α/2)
    lower = μ - z * σ / sqrt(n)
    upper = μ + z * σ / sqrt(n)
    return lower, upper
end


function waldCIHelper(df::Vector, α::Number, ind::Int)

    # this function takes the raw phenotypes data and 
    # returns the confidence intervals curves.
    CI = zeros(length(df), 2)
    for i in 1:length(df)
        data = df[i][:,ind]
        CI[i,:] .= computeWaldCI(length(data), mean(data), std(data), α)
    end
    return CI
end























