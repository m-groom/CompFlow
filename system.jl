# Functions describing the system of equations

# Load functions
include("equation_of_state.jl")

# Define the advective flux as a function of Q
function Fa(Q, γ)
    flux = zeros(3,1);
    ρ = Q[1];
    if (ρ <= 0)
        ρ = 0.0;
        u = 0.0;
        p = 0.0;
        ρE = 0.0;
    else
        u = Q[2] / Q[1];
        p = pressure(Q, γ);
        ρE = Q[3];
    end
    flux[1] = ρ * u;          # mass flux
    flux[2] = ρ * u^2 + p;    # momentum flux
    flux[3] = u * (ρE + p);   # energy flux
    return flux
end

# Define the flux Jacobian as a function of Q
function ∂F∂Q(Q, γ)
    u = Q[2]/Q[1];
    a = speedOfSound(Q, γ);
    J = [0.0 1.0 0.0;
         0.5*(γ-3.0)*u^2 (3.0-γ)*u γ-1.0;
         0.5*(γ-2.0)*u^3-a^2*u/(γ-1.0) 0.5*(3.0-2.0*γ)*u^2+a^2/(γ-1.0) γ*u];
    return J
end

# Return the eigenvalues of the flux Jacobian
function eigVal(Q, γ)
    u = Q[2]/Q[1];
    a = speedOfSound(Q, γ);
    return [u-a; u; u+a]
end

# Return the right eigenvectors of the flux Jacobian
function eigVec(Q, γ)
    ρ = Q[1]; u = Q[2]/Q[1]; ρE = Q[3];
    p = pressure(Q, γ);
    a = speedOfSound(Q, γ);
    H = (ρE + p) / ρ;
    V1 = [1.0; u-a; H-u*a];
    V2 = [1.0; u; 0.5*u^2];
    V3 = [1.0; u+a; H+u*a];
    return [V1 V2 V3]
end

# Return the left eigenvectors of the flux Jacobian
# Note: Currently assumes ideal gas EOS
function eigVecInv(Q, γ)
    u = Q[2]/Q[1];
    a = speedOfSound(Q, γ);
    V1 = [0.5*u^2+u*a/(γ-1.0) -a/(γ-1.0)-u 1.0];
    V2 = [2.0*a^2/(γ-1.0)-u^2 2.0*u -2.0];
    V3 = [0.5*u^2-u*a/(γ-1.0) a/(γ-1.0)-u 1.0];
    return 0.5 * (γ-1.0)/(a^2) * [V1; V2; V3]
end

# Convert conserved variables to primitive variables
# Note: Currently assumes ideal gas EOS
function consToPrim(Q, γ)
    ρ = Q[1]; u = Q[2]/Q[1]; ρE = Q[3];
    p = pressure(Q, γ);

    return [ρ; u; p]
end

# Convert primitive variables to conserved variables
# Note: Currently assumes ideal gas EOS
function primToCons(W, γ)
    ρ = W[1]; u = W[2]; p = W[3];
    e = internalEnergy(W, γ);
    ρE = ρ * (e + 0.5 * u^2);

    return [ρ; ρ*u; ρE]
end