function relaxSprings(x0::Vector{Float64},
                      SM::Matrix{Float64},
                      RLM::Matrix{Float64},
                      EM::Matrix{Float64};
                      gradTol::Float64 = 1e-6)::Tuple{Float64, Vector{Float64}}
    x = copy(x0)
    x = minimizeFIRE(x, SM, RLM, velocityVerlet!)    
    E = computeEnergy(x, SM, RLM, EM)
    return E, x
end

function relaxSprings(xy0::Matrix{Float64},
                      SM::Matrix{Float64},
                      RLM::Matrix{Float64},
                      EM::Matrix{Float64};
                      gradTol::Float64 = 1e-6)
    E, x = relaxSprings(xy2x(xy0),SM,RLM,EM, gradTol=gradTol)
    xy = x2xy(x)
    return E, xy
end

function relaxSprings(xy0::Matrix{Float64},
                      pheno::Vector{Matrix{Float64}};
                      gradTol::Float64 = 1e-6)
    
    E, xy = relaxSprings(xy0, pheno[1], pheno[2], pheno[3], gradTol=gradTol)
    return E, xy
end


## Relaxation algos


function minimizeFIRE(x::Vector{Float64},
                      SM::Matrix{Float64},
                      RLM::Matrix{Float64},
                      integrator!;
                      a0 = .30,
                      ashrink = .99,
                      dt0 = 0.1,
                      dt_max = 0.3,
                      dtgrow = 1.1,
                      dtshrink = 0.5,
                      delay = 2,
                      maxForce = 1e-6,
                      maxSteps = 10000)::Vector{Float64}

    a = a0
    dt = dt0
    NPG0::Int = 0
    P::Float64 = 0.0
    N::Int = 0
    converged = false

    # intialize
    F = similar(x)
    computeForce!(F,x,SM,RLM)
    v = zeros(length(x))

    while !converged
        N += 1 
        P = dot(F,v) # compute power
        if P >= 0
            NPG0 += 1
            prefac = a * norm(v) / norm(F)
            v .= (1.0 - a) .* v .+ prefac .* F
            if NPG0 > delay
                dt = min(dt*dtgrow, dt_max)
                a *= ashrink
            end
        else
            NPG0 = 0
            v .= 0
            dt *= dtshrink
            a = a0
        end
        integrator!(x, v, F, dt, SM, RLM)
        converged = isConverged(N, F, maxForce, maxSteps)
    end
    return x
end



function implicitEuler!(x, v, F, dt, SM, RLM)
    v .+= dt .* F 
    x .+= dt .* v
    computeForce!(F,x,SM,RLM)
    return nothing
end


function velocityVerlet!(x, v, F, dt, SM, RLM)
    v .+= 0.5 * dt .* F
    x .+= dt .* v
    computeForce!(F,x,SM,RLM)
    v .+= 0.5 * dt .* F
    return nothing
end

#function explicitEuler!(x, v, F, dt, SM, RLM,a)
#    x .+= dt .* v
#    computeForce!(F, x, SM, RLM)
#    v .+= dt .* F 
#    return nothing
#end; export explicitEuler!


function isConverged(stepNumber::Int,
                     force::Vector{Float64},
                     maxForce::Float64,
                     maxSteps::Int)::Bool

    stepCrit = stepNumber >= maxSteps
    forceCrit = norm(force) < maxForce
    return stepCrit || forceCrit
end



function relaxSpringsGD!(x::Vector{Float64},
                        SM::Matrix{Float64},
                        RLM::Matrix{Float64};
                        maxForce::Float64=1e-6,
                        γ::Float64=0.01)
    # return relaxed network
    g = similar(x)
    computeGradient!(g, x, SM, RLM)
    maxForceMag = norm(g)
    while maxForceMag > maxForce
        x .+= γ .* g
        computeForce!(g, x, SM, RLM)
        maxForceMag = norm(g)
    end
end


function relaxSpringsGD(x::Vector{Float64},
                        SM::Matrix{Float64},
                        RLM::Matrix{Float64};
                        maxForce::Float64=1e-6,
                        γ::Float64=0.01)::Vector{Float64}
    x = copy(x)
    relaxSpringsGD!(x,SM, RLM, maxForce=maxForce, γ=γ)
    return x
end



