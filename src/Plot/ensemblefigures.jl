  

function genConfChangeBarFig(rawData::Array{<:Array{<:Real,2},1},
                             titles::Vector{String})

    numSel = length(rawData)
    sampleSizes = zeros(numSel)
    data = Array{Array{Float64,1},1}(undef, numSel)
    for i in 1:numSel
        n = size(rawData[i],1)
        sampleSizes[i] = n
        data[i] = vec(sum(rawData[i] .== 1, dims=1) ./ n)
    end
    fig, ax = plotBarChart(data, sampleSizes, titles, "Frac. conf. changes")
    return fig
end


function genLocConfChangeBarFig(rawData, titles)

    numSel = length(rawData)
    sampleSizes = zeros(numSel)
    data = Array{Array{Float64,1},1}(undef, numSel)
    for i in 1:numSel
        n = size(rawData[i],1)
        sampleSizes[i] = n
        data[i] = vec(sum(rawData[i] .> 3, dims=1) ./ n )
    end

    fig, ax = plotBarChart(data, sampleSizes, titles, "Frac. dist. conf. changes")
    
    for a in ax; a.set_ylim([0, 0.3]);end
    return fig
end


function genAlloBarFig(rawData, titles)

    numSel = length(rawData)
    sampleSizes = zeros(numSel)
    data = Array{Array{Float64,1},1}(undef, numSel)
    for i in 1:numSel
        n = size(rawData[i],1)
        sampleSizes[i] = n
        data[i] = vec(sum(rawData[i] .> 1e-3, dims=1) ./ n )
    end
    fig, ax = plotBarChart(data, sampleSizes, titles, "Frac. Allosteric")
    return fig
end


function genAlloCCLFig(datax::Array{<:Array{<:Real,2},1},
                       datay::Array{<:Array{<:Real,2},2},
                       titles::Array{String,1};
                       labelFontSize=16)
                                                
    fig, ax = plotCorrelations(datax, datay, titles)
    numAssay, numEns = size(datay)
    for j in 1:numEns; ax[numAssay,j].set_xlabel("Conf. Change Dist.", fontsize=labelFontSize); end
    for i in 1:numAssay; ax[i,1].set_ylabel("Allostery "*titles[i], fontsize=labelFontSize); end
    for j in 1:numEns; setLimSame!(ax[:,j], "x"); end
    for a in ax
        a.set_yscale("log")
        a.set_ylim([1e-6, 1])
    end
    return fig
end




function genAlloFitFig(datax::Array{<:Array{<:Real,2},1},
                       datay::Array{<:Array{<:Real,2},2},
                       titles::Array{String,1};
                       names=["Stability", "Binding", "Specificity", "Binding", "Specificity"],
                       labelFontSize=16)    
    @assert length(names) == length(datax)
    fig, ax = plotCorrelations(datax, datay, titles)
    numAssay, numEns = size(datay)
    for j in 1:numEns; ax[numAssay,j].set_xlabel(names[j], fontsize=labelFontSize); end
    for i in 1:numAssay; ax[i,1].set_ylabel("Allostery "*titles[i], fontsize=labelFontSize); end
    for j in 1:numEns; setLimSame!(ax[:,j], "x"); end
    for a in ax
        a.set_yscale("log")
        a.set_ylim([1e-6, 1])
    end
    return fig
end



function genCCDFitFig(datax::Array{<:Array{<:Real,2},1},
                       datay::Array{<:Array{<:Real,2},1},
                       titles::Array{String,1};
                       names=["Stability", "Binding", "Specificity", "Binding", "Specificity"],
                       labelFontSize=16)
    fig, ax = plotCorrelations(datax, datay, titles)
    numEns = length(datay)
    for j in 1:numEns; ax[j].set_xlabel(names[j], fontsize=labelFontSize); end
    ax[1].set_ylabel("Conf. Change Dist.", fontsize=labelFontSize)
    setLimSame!(ax,"y")
    return fig
end


function genStabFitFig(datax::Vector{<:Matrix{<:Real}},
                       datay::Vector{<:Matrix{<:Real}},
                       titles::Vector{String};
                       names=["Stability", "Binding", "Specificity", "Binding", "Specificity"],
                       labelFontSize=16)

    @assert length(datax) == length(datay) == length(titles)
    stabMin = minimum( [ minimum(datay[i]) for i in eachindex(datay) ] ) - 0.05
    stabMax = maximum( [ maximum(datay[i]) for i in eachindex(datay) ] ) + 0.05
    fig, ax = plotCorrelations(datax, datay, titles)
    for i in 1:length(datax); ax[i].set_xlabel(names[i], fontsize=labelFontSize); end
    for i in 1:length(datax); ax[i].set_ylim([stabMin, stabMax]); end
    ax[1].set_ylabel("Stability", fontsize=labelFontSize)
    return fig     
end



#function genFitDistFig(fits::Matrix{Float64}, assay::String)
#    assert size(fits, 2) == 2
#    fig, ax = subplots(figsize(5,3))
#    formatFig!(fig)
#    mn, mx = extrema(fits)
#    ax.hist(fits[:,1], color=[.5,.5,.5])
#    ax.hist(fits[:,2], color=[.2,.7,.2], alpha=.5)
#    ax.xlabel("fitness: $assay")
#    ax.ylabel("Counts")
#    return fig
#end



function gen_q1_vs_q234Fig(ens::Vector{Network}, ens_overlaps::Vector{Vector{Float64}})

    get_ddE(Q) = allostery(Q.energies[1], Q.energies[2], Q.energies[4], Q.energies[5])
    get_q1(overlaps) = abs.(overlaps[4])
    get_q234(overlaps) = abs.(maximum(overlaps[5:7]))
    
    ddE = get_ddE.(ens)
    q1 = get_q1.(ens_overlaps)
    q234 = get_q234.(ens_overlaps)

    fig, ax = subplots(figsize=(6,4))
    formatFig!(fig)
    fig.subplots_adjust(left=0.2, bottom=0.2)
    sc = ax.scatter(q1, q234, 30, c=ddE, cmap="YlGn", edgecolors=[.8 .8 .8])
    ax.set_xlabel(L"q_1", fontsize=16)
    ax.set_ylabel(L"\mathrm{max} \{q_2,q_3, q_4\}", fontsize=16)
    cb = fig.colorbar(sc)
    cb.set_label(L"\Delta \Delta E", rotation=-90, fontsize=16, labelpad=20)
    return fig

end


