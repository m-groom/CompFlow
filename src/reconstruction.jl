# Functions for reconstructing the solution at the cell interface

# Function for performing reconstruction and returning extrapolated values
# TODO: generalise boundary conditions
function reconstruct(Q, BCs, γ, limiter = "minmod", recVars = "conserved", recType = "characteristic")
    imax = size(Q, 2);
    nVar = size(Q, 1);
    if (recVars == "primitive")
        # Convert to primitive variables
        W = zeros(nVar, imax);
        for i = 1:imax
            W[:, i] = consToPrim(Q[:, i], γ);
        end
    elseif (recVars == "conserved")
        W = Q;
    else
        error("Invalid reconstruction variables: $(recVars)")
    end
    # Piecewise linear reconstruction in space
    ΔW = zeros(nVar, imax);
    for i = 1:imax
        # Boundary conditions
        if (i==1) 
            a = zeros(3,1); # Dirichlet boundary condition
            b = W[:,i+1] - W[:,i];
        elseif (i==imax) 
            a = W[:,i] - W[:,i-1];
            b = zeros(3,1); # Dirichlet boundary condition
        else
            a = W[:,i] - W[:,i-1]; b = W[:,i+1] - W[:,i];
        end
        # Reconstruction
        if (recType == "regular")
            ΔW[:,i] = slope(a, b, limiter, recType);
        elseif (recType == "characteristic")
            # Get left and right eigenvectors
            L = eigVecInv(W[:,i], γ, recVars); R = eigVec(W[:,i], γ, recVars);
            # Project differences onto characteristic fields
            a = L * a; b = L * b;
            # Compute the slope and project back
            ΔW[:,i] = R * slope(a, b, limiter, recType);
        else
            error("Invalid reconstruction type: $(recType)")
        end
    end
    # Extrapolate to cell boundaries
    WR = W .+ 0.5 * ΔW;
    WL = W .- 0.5 * ΔW;
    if (recVars == "primitive")
        # Convert back to conserved variables
        QR = zeros(nVar, imax);
        QL = zeros(nVar, imax);
        for i = 1:imax
            QR[:, i] = primToCons(WR[:, i], γ);
            QL[:, i] = primToCons(WL[:, i], γ);
        end
    else
        QR = WR;
        QL = WL;
    end
    
    return QR, QL

end

# Compute the time derivative (Cauchy-Kovalevskaya: Q_t = -F_x) and evolve the reconstructed data
function evolve!(QR, QL, x, Δt, γ)
    imax = size(QR, 2);
    nVar = size(QR, 1);
    Qt = zeros(nVar, imax);
    for i = 1:imax
        dx = x[i+1] - x[i];
        Qt[:, i] = - (Fa(QR[:,i], γ) .- Fa(QL[:,i], γ)) / dx;
    end
    # Update the extrapolated values
    QR .= QR .+ 0.5 * Δt * Qt;
    QL .= QL .+ 0.5 * Δt * Qt;
end

# Function for calculating the slope
function slope(a, b, limiter, recType = "regular")
    if (recType == "characteristic")
        Δ = zeros(size(a));
        Δ[1] = minmod(a[1], b[1]); # Genuinely nonlinear field
        Δ[2] = superbee(a[2], b[2]); # Linearly degenerate field
        Δ[3] = minmod(a[3], b[3]); # Genuinely nonlinear field
    else
        if (limiter == "minmod")
            Δ = minmod.(a, b);
        elseif (limiter == "superbee")
            Δ = superbee.(a, b);
        elseif (limiter == "vanLeer")
            Δ = vanLeer.(a, b);
        elseif (limiter == "eno")
            Δ = eno.(a, b);
        elseif (limiter == "central")
            Δ = central.(a, b);
        elseif (limiter == "firstOrder")
            Δ = firstOrder.(a, b);
        else
            error("Invalid limiter type")
        end
    end

    return Δ
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
    β = 2.0; # Reproduces minmod for β=1 and superbee for β=2
    if (b >= 0)
        return max(0, min(β*a, b), min(a, β*b))
    else
        return min(0, max(β*a, b), max(a, β*b))
    end
end

# van Leer slope limiter
function vanLeer(a, b)
    if (abs(a) < 1e-15)
        φ = 0.0;
    else
        r = b / a;
        if (r < 0)
            φ = 0.0;
        else
            φ = 2.0 * r / (1.0 + r);
        end    
    end

    return φ * a
end

# 2nd-order ENO
function eno(a, b)
    if (abs(a) <= abs(b))
        return a
    else
        return b
    end
end

# Central difference reconstruction
function central(a, b)
    return 0.5 * (a + b)
end

# First order reconstruction
function firstOrder(a, b)
    return 0.0
end