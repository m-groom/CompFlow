# Functions for reconstructing the solution at the cell interface
# TODO: experiment with reconstruction in characteristic variables

# Load functions
include("system.jl")

# Function for performing reconstruction and returning extrapolated values
# TODO: generalise boundary conditions
function reconstruct(Q, γ)
    imax = size(Q, 2);
    nVar = size(Q, 1);
    # Convert to primitive variables
    W = zeros(nVar, imax);
    for i = 1:imax
        W[:, i] = consToPrim(Q[:, i], γ);
    end
    # Piecewise linear reconstruction in space
    ΔW = zeros(nVar, imax);
    for i = 1:imax
        if (i==1) # Boundary condition
            ΔW[:,i] = zeros(3,1); 
        elseif (i==imax) # Boundary condition
            ΔW[:,i] = zeros(3,1); 
        else 
            # TODO: generalise this to be either first or second-order
            ΔW[:,i] = minmod.(W[:,i] - W[:,i-1], W[:,i+1] - W[:,i]);
        end
    end
    # Extrapolate to cell boundaries
    WR = W .+ 0.5 * ΔW;
    WL = W .- 0.5 * ΔW;
    # Convert back to conserved variables
    QR = zeros(nVar, imax);
    QL = zeros(nVar, imax);
    for i = 1:imax
        QR[:, i] = primToCons(WR[:, i], γ);
        QL[:, i] = primToCons(WL[:, i], γ);
    end
    
    return QR, QL

end

# Compute the time derivative (Cauchy-Kovalevskaya: Q_t = -F_x) and evolve the reconstructed data
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
end

# Superbee slope limiter
function superbee(a, b)
    if (b >= 0)
        return max(0, min(2.0*a, b), min(a, 2.0*b))
    else
        return min(0, max(2.0*a, b), max(a, 2.0*b))
    end
end

# Monotonised central limiter
function mc(a, b)
    if (b >= 0)
        return max(0, min(0.5*(a+b), 2.0*a, 2.0*b))
    else
        return min(0, max(0.5*(a+b), 2.0*a, 2.0*b))
    end
end