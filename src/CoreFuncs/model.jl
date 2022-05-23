# define Energy and Force Field of model


# Energy

function computeEnergy(x::Vector{Float64},
                       SM::Matrix{Float64},
                       RLM::Matrix{Float64},
                       EM::Matrix{Float64})::Float64
    # compute the gradient of the energy function.
    # x is a vector of length 2N
    n = size(SM,1)
    energy = 0.0
    @inbounds for j in n-1:-1:1
        j2 = 2j
        xj = x[j2-1]
        yj = x[j2]
        @inbounds for i in j+1:n
            i2 = 2i
            sm = SM[i,j]
            dx = x[i2-1] - xj
            dy = x[i2] - yj
            @fastmath r = sqrt(dx^2 + dy^2)
            if sm != 0.0 || r < CUTOFF
                @fastmath energy += 0.5*sm*(r - RLM[i,j])^2
                if sm == 0.0
                    energy += 0.5*AMP*(r-CUTOFF)^2
                end
                energy += EM[i,j] 
            end
        end
    end
    return energy
end

function computeEnergy(xy::Matrix{Float64},
                       SM::Matrix{Float64},
                       RLM::Matrix{Float64},
                       EM::Matrix{Float64})::Float64
    return computeEnergy(xy2x(xy),SM, RLM, EM)
end


# Gradients and Forces


function computeGradient!(grad::Vector{Float64},
                          x::Vector{Float64},
                          SM::Matrix{Float64},
                          RLM::Matrix{Float64})
    # compute the gradient of the energy function.
    # x is a vector of length 2N
    n = size(SM,1)
    grad[2n-1] = 0.0
    grad[2n] = 0.0
    @inbounds for j in n-1:-1:1
        
        j2 = 2j
        grad[j2-1] = 0.0
        grad[j2] = 0.0
        sumx = 0.0
        sumy = 0.0
        xj = x[j2-1]
        yj = x[j2]
        @inbounds for i in j+1:n
            sm = SM[i,j]
            i2 = 2i
            dx = xj - x[i2-1]
            dy = yj - x[i2]
            @fastmath r = sqrt(dx^2 + dy^2)
            if sm != 0.0 || r < CUTOFF
                fp = sm * (r - RLM[i,j])/r 
                if r < CUTOFF
                    fp += AMP * (r - CUTOFF)/r
                end
                fx = fp*dx
                fy = fp*dy
                sumx += fx
                sumy += fy
                grad[i2-1] -= fx
                grad[i2] -= fy
            end
        end
        grad[j2-1] = sumx
        grad[j2] = sumy
    end
    return nothing
end 


function computeForce!(F::Vector{Float64},
                       x::Vector{Float64},
                       SM::Matrix{Float64},
                       RLM::Matrix{Float64})
    computeGradient!(F, x, SM, RLM)
    F .*= -1 # convert the gradient to the force
    return nothing
end


function computeForce(x::Vector{Float64},
                      SM::Matrix{Float64},
                      RLM::Matrix{Float64})
    F_out = similar(x)
    computeForce!(F_out, x, SM, RLM)
    return F_out
end


function computeForce(xy::Matrix{Float64},
                      SM::Matrix{Float64},
                      RLM::Matrix{Float64})
    force = computeForce(xy[:], SM, RLM)
    force = reshape(force, size(xy))
    return force
end


##### Hessian ##############



function computeHessian!(H::Matrix{Float64},
                         x::Vector{Float64},
                         SM::Matrix{Float64},
                         RLM::Matrix{Float64})
    # compute Hessian matrix of Hamiltonian
    # returns hessian
    N = size(SM,1)
    H .= 0
    for i in 1:N
        sxx = 0.0
        syy = 0.0 
        sxy = 0.0
        i2 = 2i
        xi = x[i2-1]
        yi = x[i2]
        for j in 1:N
            i == j && continue
            j2 = 2j
            dx = xi - x[j2-1]
            dy = yi - x[j2]
            sm = SM[i,j]
            r = sqrt(dx^2+dy^2)
            if sm != 0.0 || r < CUTOFF
                l = RLM[i,j]
                Exx = sm*( dx^2*(r-l)/r^3 - dx^2/r^2 - (r-l)/r )
                Eyy = sm*( dy^2*(r-l)/r^3 - dy^2/r^2 - (r-l)/r )
                Exy = sm*( dx*dy*(r-l)/r^3 - dx*dy/r^2 )
                if r < CUTOFF
                    Exx += AMP*( dx^2*(r-CUTOFF)/r^3 - dx^2/r^2 - (r-CUTOFF)/r )
                    Eyy += AMP*( dy^2*(r-CUTOFF)/r^3 - dy^2/r^2 - (r-CUTOFF)/r )
                    Exy += AMP*( dx*dy*(r-CUTOFF)/r^3 - dx*dy/r^2 )
                end   
                sxx -= Exx
                syy -= Eyy
                sxy -= Exy
                H[2i-1,2j-1] = Exx
                H[2i,2j] = Eyy
                H[2i,2j-1] = Exy
                H[2i-1,2j] = Exy
            end
        end
        H[2i-1,2i-1] = sxx
        H[2i,2i] = syy
        H[2i,2i-1] = sxy
        H[2i-1,2i] = sxy
    end
    return nothing
end

old2new(x) = reshape(x,:,2)'[:]
new2old(x) = reshape(x,2,:)'[:]


#function compHess!(H, x, SM, RLM)
#    computeHessian!(H, new2old(x), SM, RLM)
#end; export compHess

function computeHessian(x::Vector{Float64}, SM::Matrix{Float64}, RLM::Matrix{Float64})
    H = zeros(2*size(SM,1), 2*size(SM,1))
    computeHessian!(H, x, SM, RLM)
    return H
end

function computeHessian(Q::Network)
    return computeHessian(Q.structs[1][:], Q.pheno[1], Q.pheno[2])
end


