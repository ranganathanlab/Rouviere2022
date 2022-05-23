

function evolve(Q0::Network,
                assay::String,
                evoT::Float64,
                numEvoSteps::Int;
                GS_settings::Tuple{Int, Int}=(10,100),
                E_thresh::Number=0.0001,
                streekThresh::Int=25,
                α::Int=5,# mutiplied by numTrials for extra confidence
                showProgress::Bool=false)
    # internal parameters
    Q = deepcopy(Q0)
    Q_old = deepcopy(Q)
    Qs = Network[]; 
    energiesList = Matrix{Float64}[]
    fits = Float64[]
    evoSteps = Int[]
    ligDex = setLigDex(assay)
    special_GS_settings = (α*GS_settings[1], GS_settings[2])
    # Initial state
    findGroundStates!(Q, setLigDex(assay), special_GS_settings)
    fitness = computeFitness(Q, assay)
    oldFitness = fitness
    saveState!(energiesList, fits, evoSteps, Q.energies, fitness, 0)
    push!(Qs, deepcopy(Q))

    i = 0; streek = 0
    while (i < numEvoSteps) || (streek < streekThresh)
        i += 1; streek += 1
        showProgress && displayProgress(numEvoSteps, i) # display progress
        # finish by locally optimizing
        (i > numEvoSteps - streekThresh) && (evoT = zero(evoT))

        # check to make sure you are in ground state
        findGroundStates!(Q, ligDex, GS_settings)
        fitness_test = computeFitness(Q, assay)
        if fitness_test < fitness - E_thresh
            i -= 1
            Q = deepcopy(Q_old)
            findGroundStates!(Q, ligDex, GS_settings)
            fitness = computeFitness(Q, assay)
            if length(fits) > 1
                pop!(fits); pop!(energiesList); pop!(evoSteps)
                fits[end] = fitness
                energiesList[end] = copy(Q.energies)
            else
                fits[1] = fitness; energiesList[1] = copy(Q.energies) 
            end
        end
        # mutate and compute new fitness
        Q_wildtype = deepcopy(Q)
        randomMutateNode!(Q)
        findGroundStates!(Q, ligDex, GS_settings)
        newFitness = computeFitness(Q, assay)
        # accept mutant with monte carlo condition
        if acceptMC(fitness - newFitness, evoT)
            oldFitness = fitness
            fitness = newFitness
            Q_old = deepcopy(Q_wildtype)
            saveState!(energiesList, fits, evoSteps, Q.energies, fitness, i)
        else
            Q = deepcopy(Q_wildtype)
        end

        # if you think you are done, double check
        if  (i >= numEvoSteps) && (streek >= streekThresh)
            energies, structs = findGroundStates(Q, ligDex, special_GS_settings)
            if any((energies .- Q.energies[ligDex]) .< -E_thresh) # if you were wrong try again
                println("you where not in the GS")
                streek = 0
                Q.energies[ligDex] .= energies
                Q.structs[ligDex] .= deepcopy.(structs)
                fitness = computeFitness(Q, assay)
            else # if you were right, save and end
                push!(Qs,deepcopy(Q))
            end
        end
    end # exit evolution loop
    findGroundStates!(Qs[end], ligDex,special_GS_settings)
    energiesList[end] .= Qs[end].energies
    return Qs, energiesList, fits,  evoSteps
end




function saveState!(energiesList, fits, evoSteps, energies, fitness, evoStep )
    push!(energiesList, copy(energies))
    push!(fits, fitness)
    push!(evoSteps, evoStep)
end


function displayProgress(totalSteps::Int, currentStep::Int)
    print("$(computeProgress(currentStep, totalSteps))%\r")
end


function setLigs!(Q::Network, ligsList, ligIndex::Int)
    ligIndex = (ligIndex % length(ligsList)) + 1
    Q.ligs = ligsList[ligIndex]
    return ligIndex
end

function setLigs!(Qs::Vector{Network}, ligsList, ligIndex::Int)
    ligIndex = (ligIndex % length(ligsList)) + 1
    for Q in Qs
        Q.ligs = ligsList[ligIndex]    
    end
    return ligIndex
end

###########################################

function run1(S::Dict, seed::UInt)
    # Generate Random Network and evolve it.

    Random.seed!(seed) # sets the seed for the type table nothing else.
    Q = buildNetworkFromSettings(S)

    Qs, energies, fits, evoSteps = evolve(Q, S["assay"], S["evoT"], S["numEvoSteps"],
                                          GS_settings = S["GS_settings"])
    return Qs, energies, fits, evoSteps
end


function run1(S::Dict, Q::Network)
    # Evolve the network Q 
    Qs, energies, fits, evoSteps = evolve(Q, S["assay"], S["evoT"], S["numEvoSteps"],
                                          GS_settings = S["GS_settings"])
    
    return Qs, energies, fits, evoSteps
end




function runMany(S::Dict, numReplicates::Int, seeds::Vector{UInt})
    # Generate numReplicates random networks and evolve them.
    out = pmap(x -> run1(S, x), seeds)
    ensemble = Array{Network,2}(undef, numReplicates, 2)
    ensemble_energies = Vector{Vector{Matrix{Float64}}}(undef, numReplicates)
    ensemble_fits = Vector{Vector{Float64}}(undef, numReplicates)
    ensemble_evoSteps = Vector{Vector{Int}}(undef, numReplicates)
    for i in 1:numReplicates
        Qs, energies, fits, evoSteps = out[i]
        ensemble[i,:] = Qs
        ensemble_energies[i] = energies
        ensemble_fits[i] = fits
        ensemble_evoSteps[i] = evoSteps
    end
    return ensemble, ensemble_energies, ensemble_fits, ensemble_evoSteps
end

function runMany(S::Dict,
                 ensemble::Array{Network})
    # evolve the networks provided in ensemble.
    numReplicates = length(ensemble)
    out = pmap(x -> run1(S, x), ensemble)
    ensemble = Array{Network,2}(undef, numReplicates, 2)
    ensemble_energies = Vector{Vector{Matrix{Float64}}}(undef, numReplicates)
    ensemble_fits = Vector{Vector{Float64}}(undef, numReplicates)
    ensemble_evoSteps = Vector{Vector{Int}}(undef, numReplicates)
    for i in 1:numReplicates
        Qs, energies, fits, evoSteps = out[i]
        ensemble[i,:] = Qs
        ensemble_energies[i] = energies
        ensemble_fits[i] = fits
        ensemble_evoSteps[i] = evoSteps
    end
    return ensemble, ensemble_energies, ensemble_fits, ensemble_evoSteps
end










################### maybe move this to families #################

function evolveFamily(Q0::Network,
                      assay::String,
                      evoT::Float64,
                      numEvoSteps::Int,
                      numTrials::Int,
                      famSize::Int)
    
    function kernel(Q)
        Q = deepcopy(Q) # detach from outside 
        Qs, fits, evoSteps = evolveLite(Q, assay, evoT, numEvoSteps, numTrials=numTrials, saveIntermediates=true)
        Ls = makeLilNetwork.(Qs)
        return Ls, fits, evoSteps
    end
    return pmap(x -> kernel(Q0), 1:famSize)
end


function formatFamily(evolveFamilyOutput, numEvoSteps)
    # reformatts to have 
    famSize = length(evolveFamilyOutput)
    fitnesses = zeros(numEvoSteps+1, famSize)
    sequences = Array{Vector{Int}, 2}(undef, numEvoSteps+1, famSize)
    seqs_temp = [[i.geno for i in j[1]] for j in evolveFamilyOutput] 
    for (i, tri) in enumerate(evolveFamilyOutput)
        Ls, fits, evoSteps = tri
        listOseqs = seqs_temp[i]
        sequences[:,i] = makeEvoTrace(evoSteps, listOseqs, numEvoSteps)
        fitnesses[:,i] = makeEvoTrace(evoSteps, fits, numEvoSteps)
    end
    return fitnesses, sequences    
end


function makeEvoTrace(evoSteps::Vector{Int},
                      things::Vector{<:Number},
                      numEvoSteps::Int)
    #@assert length(evoSteps) == length(things)
    numEvoSteps=Int(numEvoSteps)
    evoTrace = Vector{eltype(things)}(undef, numEvoSteps+1)
    for i in 1:length(evoSteps)-1
        step1 = evoSteps[i]
        step2 = evoSteps[i+1]
        for j in step1+1:step2
            evoTrace[j] = things[i]
        end
    end
    
    for j in evoSteps[end]+1:numEvoSteps+1
        evoTrace[j] = things[end]
    end
    return evoTrace
end



## population evolution ###


function evolvePop(Qs::Vector{Network},
                   assay::String,
                   p::Number,
                   numEvoSteps::Int,
                   τ::Int,
                   ligsList::Vector{Vector{Ligand}};
                   GS_settings::Tuple{Int, Int}=(10,100),
                   )

    popSize = length(Qs)
    num2kill = Int(round(length(Qs)*p))
    ligDex = setLigDex(assay)
    fitsArray = zeros(popSize, numEvoSteps)

    Qs = pmap(findGroundStates2, Qs, fill(ligDex, popSize), fill(GS_settings, popSize))
    fits = map(x -> computeFitness(x, assay), Qs)
    fits0 = copy(fits)
    fitsArray[:,1] = fits 

    ligIndex = setLigs!(Qs, ligsList, length(ligsList))

    for i in 1:numEvoSteps

        losers = pickLosers(fits, num2kill)
        winners = pickWinners(losers, popSize)
        Qs[losers] = deepcopy.(Qs[winners])
        randomMutateNode!.(Qs[losers])

        # swap the ligands
        if ((i-1) % τ == 0) && !(τ == numEvoSteps)
            ligIndex = setLigs!(Qs, ligsList, ligIndex)
        end 

        Qs = pmap( findGroundStates2, Qs, fill(ligDex, popSize), fill(GS_settings, popSize))
        fits = map(x->computeFitness(x, assay), Qs)
        fitsArray[:,i] = fits 
    end
    return Qs, fitsArray
end


function pickLosers(fits::Vector, num2kill::Int)
    return sortperm(fits)[1:num2kill]
end


function pickWinners(losers::Vector, popSize::Int)
    return rand(setdiff(1:popSize,losers), length(losers))
end

function evolvePopWrapper(S::Dict)
    Qs0 = [buildNetworkFromSettings(S) for i in 1:S["popSize"]]
    Qs1, fitsArray = evolvePop(Qs0, S["assay"], S["p"], S["numEvoSteps"],
                              S["τ"],  S["ligsList"],
                              GS_settings=S["GS_settings"]);

    Qs0 = pmap( findGroundStates2, Qs0, fill(setLigDex("Specificity"), S["popSize"]), fill(S["GS_settings"], S["popSize"]));
    Qs1 = pmap( findGroundStates2, Qs1, fill([CartesianIndex(3, 1)], S["popSize"]), fill(S["GS_settings"], S["popSize"]));
    return [Qs0 Qs1], fitsArray
end











########################################
### Flux ###############################



function evolveFlux(Q0::Network,
                    assay::String,
                    evoT::Float64,
                    numEvoSteps::Int,
                    τ::Int,
                    ligsList::Vector{Vector{Ligand}};
                    GS_settings::Tuple{Int, Int}=(10,100),
                    E_thresh::Number=0.0001,
                    streekThresh::Int=25,
                    α::Int=5,# mutiplied by numTrials for extra confidence
                    showProgress::Bool=false,
                    )
    
    # internal parameters
    Q = deepcopy(Q0)
    Q_old = deepcopy(Q)
    Qs = Network[]; 
    energiesList = Matrix{Float64}[]
    fits = Float64[]
    evoSteps = Int[]
    ligDex = setLigDex(assay)
    special_GS_settings = (α*GS_settings[1], GS_settings[2])

    # Initial state
    findGroundStates!(Q, setLigDex("All"), special_GS_settings)
    fitness = computeFitness(Q, assay)
    oldFitness = fitness
    saveState!(energiesList, fits, evoSteps, Q.energies, fitness, 0)
    push!(Qs, deepcopy(Q))
    ligIndex = setLigs!(Q, ligsList, length(ligsList))

    i = 0
    while i < numEvoSteps
        i += 1
        showProgress && displayProgress(numEvoSteps, i) # display progress
        
        # check to make sure you are in ground state
        findGroundStates!(Q, ligDex, GS_settings)
        fitness_test = computeFitness(Q, assay)
        if fitness_test < fitness - E_thresh
            i -= 1
            Q = deepcopy(Q_old)
            findGroundStates!(Q, ligDex, GS_settings)
            fitness = computeFitness(Q, assay)
           
            if length(fits) > 1
                pop!(fits); pop!(energiesList); pop!(evoSteps)
                fits[end] = fitness
                energiesList[end] = copy(Q.energies)
            else
                fits[1] = fitness; energiesList[1] = copy(Q.energies) 
            end
        end

        # swap the ligands
        if ((i-1) % τ == 0) && !(τ == numEvoSteps)
            ligIndex = setLigs!(Q, ligsList, ligIndex)
            findGroundStates!(Q, ligDex, GS_settings)
            fitness = computeFitness(Q, assay)
        end 

        # mutate and compute new fitness
        Q_wildtype = deepcopy(Q)
        randomMutateNode!(Q)
        findGroundStates!(Q, ligDex, GS_settings)
        newFitness = computeFitness(Q, assay)
        
        # accept mutant with monte carlo condition
        if acceptMC(fitness - newFitness, evoT)
            oldFitness = fitness
            fitness = newFitness
            Q_old = deepcopy(Q_wildtype)
            saveState!(energiesList, fits, evoSteps, Q.energies, fitness, i)
        else
            Q = deepcopy(Q_wildtype)
        end
    end # exit evolution loop
    push!(Qs, deepcopy(Q))
    findGroundStates!(Qs[end], setdiff(setLigDex("All"), ligDex), special_GS_settings)
    energiesList[end] .= Qs[end].energies
    return Qs, energiesList, fits,  evoSteps
end

function runManyFlux(S::Dict, numReplicates::Int, seeds::Vector{UInt})
    out = pmap(x -> run1Flux(S, x), seeds)
    ensemble = Array{Network,2}(undef, numReplicates, 2)
    ensemble_energies = Vector{Vector{Matrix{Float64}}}(undef, numReplicates)
    ensemble_fits = Vector{Vector{Float64}}(undef, numReplicates)
    ensemble_evoSteps = Vector{Vector{Int}}(undef, numReplicates)
    for i in 1:numReplicates
        Qs, energies, fits, evoSteps = out[i]
        ensemble[i,:] = Qs
        ensemble_energies[i] = energies
        ensemble_fits[i] = fits
        ensemble_evoSteps[i] = evoSteps
    end
    return ensemble, ensemble_energies, ensemble_fits, ensemble_evoSteps
end



function run1Flux(S::Dict, seed::UInt)
    Random.seed!(seed) # sets the seed for the type table nothing else.
    Q = buildNetworkFromSettings(S)
    Qs, energies, fits, evoSteps = evolveFlux(Q, S["assay"],
                                  S["evoT"], S["numEvoSteps"],
                                  S["τ"],  S["ligsList"],
                                  GS_settings=S["GS_settings"])
    return Qs, energies, fits, evoSteps
end
