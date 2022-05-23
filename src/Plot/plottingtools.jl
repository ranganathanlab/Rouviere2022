

function clr(colorChoice)
    # colors are going to be used with scatter() then
    # they must be 2-d arrays
    colors = Dict()
    colors["blue"] = [.2, .3, 1]
    colors["green"] = [.2, .7, .2]
    colors["red"] = [.8, .1, .1]
    colors["black"] = [0, 0, 0]
    colors["bblack"] = [0 0 0]
    colors["grey"] = [.7, .7, .7] #
    colors["ggrey"] = [.7 .7 .7]
    colors["darkgrey"] = [.4 .4 .4]
    colors["b1"] = [.4 .4 .9 0.5]
    colors["g1"] = [.2 .6 .2 0.5]
    return colors[colorChoice]
end



function formatFig!(fig)
    # this function sets my plotting preferences
    ax = fig.axes
    for i in 1:length(ax)
        ax[i].set_title("",fontsize=16)
        ax[i].set_xlabel("",fontsize=14)
        ax[i].set_ylabel("",fontsize=14)
        for axs in ["top","bottom","right","left"]
            ax[i].spines[axs].set_linewidth(1.3)
        end
        ax[i].xaxis.set_tick_params(width=1.3)
        ax[i].yaxis.set_tick_params(width=1.3)
    end
    return nothing
end


#function jet(numColors)
#    # uses PyCall
#    maps = pyimport("matplotlib.cm")
#    rawColors = maps.jet(0:255)[:,1:3]
#    r =Int.(round.(range(1, stop=256 , length = numColors)))
#    colors = rawColors[r,:]
#    return colors
#end; export jet


function Reds(numColors)
    colorArray = ones(numColors,3)
    grad = linspace(1,0,numColors)
    colorArray[:,2:3] .= [grad grad]
    return colorArray
end



function makeColorbar!(fig; cm="Reds", low=0, high=1, x=0, y=0, wid=.02, len=.5)
    clrmp = get_cmap(cm)
    mpl = pyimport("matplotlib")
    norm = mpl.colors.Normalize(vmin=low,vmax=high)
    sm = mpl.cm.ScalarMappable(cmap=clrmp, norm=norm)
    sm.set_array([])
    cbaxes = fig.add_axes([x, y, wid, len])
    fig.colorbar(sm, ticks=[low,high], cax=cbaxes)
    return nothing
end



function makeAllNetworkAxSame!(axArray)
    # designed for plotting elastic network structures
    xmin=Inf
    xmax=-Inf
    ymin=Inf
    ymax=-Inf
    padPercent = .01

    for ax in axArray
        low, high = ax.get_xlim()
        if low < xmin; xmin=low;end
        if high > xmax; xmax=high;end

        low, high = ax.get_ylim()
        if low < ymin; ymin=low;end
        if high > ymax; ymax=high;end
    end
    
    padx = abs(xmax-xmin)*padPercent
    pady = abs(ymax-ymin)*padPercent

    for ax in axArray
        ax.set_xlim( [xmin-padx, xmax+padx] )
        ax.set_ylim( [ymin-pady, ymax+pady] )
        ax.set_aspect("equal",adjustable="box")
        ax.axis("off")
    end    
    return nothing
end



function setLimSame!(axArray, axis)
    # sets limits of axes to the min and max for each ax object

    if axis!="x" && axis!="y"
        error("choose an axis, either \"x\" or \"y\"")
    end

    r = []
    for a in axArray
        if axis == "x"
            append!(r, a.get_xlim())
        elseif axis == "y"
            append!(r, a.get_ylim())
        end
    end
    r = extrema(r)
    for a in axArray
        if axis == "x"
            a.set_xlim(r)
        elseif axis == "y"
            a.set_ylim(r)
        end
    end
    return nothing
end



function genBlankFig(txt)
    fig, ax = subplots(1,1,figsize=(10,2))
    ax.axis("off")
    text(0,.5, txt, fontsize=50)
    return fig
end


function jitter!(ax,
                x::Number,
                data::Vector{<:Number},
                jit::Number;
                ps::Number=5,
                c="k",
                alpha::Number=1)

    xs = jit .* (rand(length(data)) .- 0.5) .+ x
    ax.scatter( xs, data, ps, c=c, alpha=alpha )
end


function jitter!(ax,
                 x::AbstractVector{<:Number},
                 data::Vector{<:Vector{<:Number}},
                 jit::Number;
                 ps::Number=5,
                 c="k",
                 alpha::Number=1)

    @assert length(x) == length(data)

    for i in 1:length(data)
        slice=data[i]
        xx = x[i]
        jitter!(ax, xx, slice, jit, ps=ps, c=c, alpha=alpha)
    end
end
                 




























