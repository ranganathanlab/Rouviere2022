

function computeEvoTraces(Q::Network,
                          assay::String;
                          evoT::Float64=0.002,
                          numEvoSteps::Int=25,
                          numTrials::Int=1000,
                          numReps::Int=25)::Matrix{Float64}
    # evo traces are stored in cols of returned matrix 
    function reevolve()
        Qs, fits, evoSteps = evolveLite(Q, assay, evoT,
            numEvoSteps, numTrials=numTrials, saveIntermediates=true)
        evoTrace = makeEvoTrace(evoSteps, fits, numEvoSteps)
        return evoTrace
    end
    out = map(x -> reevolve(), 1:numReps)
    return hcat(out...) 
end


#function makeEvoTrace(evoSteps::Vector{Int},
#                      things::Vector{T},
#                      numEvoSteps::Int) where T
#    @assert length(evoSteps) == length(things)
#    evoTrace = Vector{eltype(things)}(udef, numEvoSteps+1)
#    for i in 1:length(evoSteps)-1
#        step1 = evoSteps[i]
#        step2 = evoSteps[i+1]
#        evoTrace[step1+1:step2] = things[i]
#    end
#    evoTrace[evoSteps[end]+1:end] = things[end]
#    return evoTrace
#end




function ensembleComputeEvoTraces(ensemble::Matrix{Network},
                                  assay::String;
                                  evoT::Float64=0.002,
                                  numEvoSteps::Int=25,
                                  numReps::Int=25)

    return pmap( x -> computeEvoTraces(x, assay, evoT=evoT, numEvoSteps=numEvoSteps, numReps=numReps), ensemble)
end


function sweepComputeEvoTraces(BigData::Array{Any},
                               assay::String,
                               evoT::Float64,
                               numEvoSteps::Int=25,
                               numReps::Int=25)

    return map( x -> computeEvoTraces(x, assay, evoT, numEvoSteps=numEvoSteps, numReps=numReps), BigData)
end




function ensembleComputeEvol(ensembleEvoTraces::Array{Any,2})
    ensembleEvol = computeEvol.(ensembleEvoTraces)
    return ensembleEvol
end



#function computeEvol(evoTraces::Array{Array{Float64,1},1})
#    avgTrace = mean( hcat(evoTraces...) , dims=2)
#    return avgTrace[end] - avgTrace[1]
#end










