# individual make all figs


function makeAll_energies(ensemble,
                          dir::String)

    function makeNsave(Q, num, dir)
        fig, name = genEnergyTickFig(Q)
        savve(fig, name, dir, num)
    end

    N = size(ensemble,1)
    numbers = 1:N
    map(makeNsave, ensemble, numbers, fill(dir,N))
end



function makeAll_energyVsConfCoor(ensemble,
                                    dir::String)

    function makeNsave(Q, num, dir)
        fig, name = genEnergyVsConfCoorFig(Q)
        savve(fig, name, dir, num)
    end

    N = size(ensemble,1)
    numbers = 1:N
    map(makeNsave, ensemble, numbers, fill(dir,N))
end


function makeAll_structs(ensemble,
                         dir::String)

    function makeNsave(Q, num, dir)
        fig, name = genStructsFig(Q)
        savve(fig, name, dir, num)
    end
    N = size(ensemble,1)
    numbers = 1:N
    map(makeNsave, ensemble, numbers,fill(dir,N))
end


function makeAll_structsAllo(ensemble,
                             dir::String)

    function makeNsave(Q, num, dir)
        fig, name = genStructsAlloFig(Q)
        savve(fig, name, dir, num)
    end
    N = size(ensemble,1)
    numbers = 1:N
    map(makeNsave, ensemble, numbers,fill(dir,N))
end


function makeAll_strain(ensemble,
                         dir::String)

    function makeNsave(Q, num, dir)
        fig, name = genStrainMapFig(Q)
        savve(fig, name, dir, num)
    end
    N = size(ensemble,1)
    numbers = 1:N
    map(makeNsave, ensemble, numbers, fill(dir,N))
end


function makeAll_strainAllo(ensemble,
                         dir::String)
    function makeNsave(Q, num, dir)
        fig, name = genStrainAlloMapFig(Q)
        savve(fig, name, dir, num)
    end

    N = size(ensemble,1)
    numbers = 1:N
    map(makeNsave, ensemble, numbers, fill(dir,N))
end


function makeAll_mutSensMap(ensemble,
                            df,
                            dir::String,
                            assay::String)

    function makeNsave(Q, df, num, assay, dir)
        fig, name = genMutSensMapAbsFig(Q, df, assay)
        savve(fig, name, dir, num)
    end
    N = size(ensemble,1)
    numbers = 1:N
    map(makeNsave, ensemble, df, numbers, fill(assay,N), fill(dir,N))
end




function makeAll_DMS(df,
                     dir::String,
                     assay::String)

    function makeNsave(df, num, assay, dir)
        fig, name = genDMSFig(df, assay)
        savve(fig, name, dir, num)
    end

    N = size(df,1)
    numbers = 1:N
    map((x,y) -> makeNsave( x, y, assay, dir),  df, numbers)
end


function makeAll_mutSensDist(df,
                             dir::String,
                             assay::String)

    function makeNsave(df, num, assay, dir)
        fig, name = genMutSensDistFig(df, assay)
        savve(fig, name, dir, num)
    end

    N = size(df,1)
    numbers = 1:N
    map(makeNsave,  df, numbers, fill(assay,N), fill(dir,N))
end


function makeAll_mutSensRDF(ensemble,
                            df,
                            dir::String,
                            assay::String)

    function makeNsave(Q, df, num, assay, dir)
        ASNs = [Q.sites[1][1][1], Q.sites[1][1][2]]
        fig, name = genMutSensRDFFig(df, Q.A, ASNs, assay);
        savve(fig, name, dir, num)
    end

    N = size(df,1)
    numbers = 1:N
    map(makeNsave, ensemble, df, numbers, fill(assay,N), fill(dir,N))
end


function makeAll_mutSens(ensemble,
                         GSE_WT,
                         GSE_DMS,
                         dir,
                         assay)

    # make a mut sens figs for particular assay
    params  = constructParams(ensemble[1], assay)
    df = ensembleComputeFitSens(GSE_WT, GSE_DMS, params)
    
    makeAll_mutSensMap(ensemble, df, dir, assay)
    #makeAll_DMS(df, dir, assay)
    makeAll_mutSensDist(df, dir, assay)
    #makeAll_mutSensRDF(ensemble, df, dir, assay)
    return nothing
end



function makeAll_StabBind(ensemble, GSE_WT, GSE_DMS, dir)

    function makeNsave(Q, GSE_WT, GSE_DMS, num, dir)
        params = constructParams(Q, "Stability")
        stab = computeSingleMutsFit(GSE_DMS, params)
        stabWT = computeFitness(GSE_WT, params)

        params = constructParams(Q, "Binding")
        bind = computeSingleMutsFit(GSE_DMS, params)
        bindWT = computeFitness(GSE_WT, params)

        # get rid of neutral mutations.
        s = Float64[]
        b = Float64[]
        for i in 1:size(stab,1), j in 1:size(stab,2)
            if j != Q.geno[i]
                push!(s, stab[i,j])
                push!(b, bind[i,j])
            end
        end
                
        fig, name = genStabBindScatterHistFig(s,b,stabWT, bindWT, unfoldedEnergy=Q.unfoldedEnergy);
        savve(fig, name, dir, num)
    end

    N = size(ensemble,1)
    numbers = 1:N
    map(makeNsave, ensemble, GSE_WT, GSE_DMS, numbers, fill(dir,N))

end



function makeAll_ligScape(ensemble,
                          ligScape_energies,
                          ligScape_structs,
                          ligLengths,
                          dir::String)

    function makeNsave(Q, energies, structs, ligLengths, num, dir)
        fig, name = genLigScapeFig(energies, structs, ligLengths, Q);
        savve(fig, name, dir, num)
    end
    numbers = 1:length(ensemble)
    map( (v1,v2,v3, v4) -> makeNsave(v1, v2, v3, ligLengths, v4, dir), ensemble, ligScape_energies, ligScape_structs, numbers)
end


function makeAll_strainVsSoftMode(ensemble,
                              dir::String)

    function makeNsave(Q, num, dir)
        fig, name, fig2, name2 = genStrainVsSoftModeFig(Q)
        savve(fig, name, dir, num)
        savve(fig2, name2, dir, num)
    end

    N = size(ensemble,1)
    numbers = 1:N
    map(makeNsave, ensemble, numbers, fill(dir,N))
end

#########################################

function savve(fig, name, dir, number)
    # save figure in proper directory with proper
    # name

    num = @sprintf("%03.0f", number)        
    path =  dir*num*"_"*name*".svg"
    fig.savefig(path, format="svg", bbox_inches = "tight")
    close(fig)
    return nothing
end

