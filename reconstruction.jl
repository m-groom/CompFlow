# Functions for reconstructing the solution at the cell interface
include("system.jl")

# Function for performing reconstruction and returning
# TODO: generalise boundary conditions
function reconstruct(Q)
    imax = size(Q, 2);
    nVar = size(Q, 1);
    ΔQ = zeros(nVar, imax);
    QR = zeros(nVar, imax);
    QL = zeros(nVar, imax);
    for i = 1:imax
        if (i==1) # Boundary condition
            ΔQ[:,i] = zeros(3,1); 
        elseif (i==imax) # Boundary condition
            ΔQ[:,i] = zeros(3,1); 
        else
            ΔQ[:,i] = minmod.(Q[:,i] - Q[:,i-1], Q[:,i+1] - Q[:,i]);
        end
    end
    # Piecewise linear reconstruction in space
    # TODO: generalise this to be either first or second-order
    QR = Q .+ 0.5 * ΔQ;
    QL = Q .- 0.5 * ΔQ;
    
    return QR, QL

end

# Compute the time derivative (Cauchy-Kovalevskaya): Q_t = -F_x
function evolve(QR, QL, x, dt, γ)
    imax = size(QR, 2);
    nVar = size(QR, 1);
    Qt = zeros(nVar, imax);
    for i = 1:imax
        dx = x[i+1] - x[i];
        Qt[:, i] = - (Fa(QR[:,i], γ) .- Fa(QL[:,i], γ)) / dx;
    end
    # Update the extrapolated values
    QR = QR .+ 0.5 * dt * Qt;
    QL = QL .+ 0.5 * dt * Qt;

    return QR, QL

end

# Minmod slope limiter
function minmod(a, b)
    if (b >= 0)
        return max(0, min(a, b))
    else
        return min(0, max(a, b))
    end
    # return max(0, min(a, b))
end

# Superbee slope limiter
function superbee(a, b)
    if (b >= 0)
        return max(0, min(2.0*a, b), min(a, 2.0*b))
    else
        return min(0, max(2.0*a, b), max(a, 2.0*b))
    end
    # return max(0, min(2.0*a, b), min(a, 2.0*b))
end