# Functions describing the system of equations
include("equation_of_state.jl")

# Define the (inviscid) flux as a function of Q
function Fx(Q, γ)
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
function ∂Fx(Q, γ)
    u = Q[2]/Q[1];
    a = speedOfSound(Q, γ);
    J = [0.0 1.0 0.0;
         0.5*(γ-3.0)*u^2 (3.0-γ)*u γ-1.0;
         0.5*(γ-2.0)*u^3-a^2*u/(γ-1.0) 0.5*(3.0-2.0*γ)*u^2+a^2/(γ-1.0) γ*u;];
    return J
end

# Return the eigenvalues of the flux Jacobian
function eigval(Q, γ)
    u = Q[2]/Q[1];
    a = speedOfSound(Q, γ);
    return [u-a; u; u+a]
end

# Return the eigenvectors of the flux Jacobian
function eigvec(Q, γ)
    ρ = Q[1]; u = Q[2]/Q[1]; ρE = Q[3];
    p = pressure(Q, γ);
    a = speedOfSound(Q, γ);
    H = (ρE + p) / ρ;
    V1 = [1.0; u-a; H-u*a];
    V2 = [1.0; u; 0.5*u^2];
    V3 = [1.0; u+a; H+u*a];
    return [V1 V2 V3]
end

# Convert conserved variables to primitive variables
function primitive(Q, γ)
    ρ = Q[1]; u = Q[2]/Q[1]; ρE = Q[3];
    p = pressure(Q, γ);

    return [ρ; u; p]
end

# Convert primitive variables to conserved variables
function conserved(W, γ)
    ρ = W[1]; u = W[2]; p = W[3];
    e = internalEnergy(W, γ);
    ρE = ρ * (e + 0.5 * u^2);

    return [ρ; ρ*u; ρE]
end