# Functions for computing the solution of the Riemann problem
include("equation_of_state.jl")
include("newton.jl")
include("system.jl")

# Compute the exact solution of the Riemann problem (from p152 of Toro)
# TODO: replace primitve and conserved variable calculations with functions
function exactRiemannSolver(QL, QR, ξ, γ)
    # Calculate primitive variables
    DL = QL[1];        # left density
    DR = QR[1];        # right density
    UL = QL[2]/QL[1];  # left velocity
    UR = QR[2]/QR[1];  # right velocity
    PL = pressure(QL, γ);   # left pressure
    PR = pressure(QR, γ);   # right pressure

    # Compute γ related constants
    G1 = (γ - 1.0) / (2.0 * γ);
    G2 = (γ + 1.0) / (2.0 * γ);
    G3 = 2.0 * γ / (γ - 1.0);
    G4 = 2.0 / (γ - 1.0);
    G5 = 2.0 / (γ + 1.0);
    G6 = (γ - 1.0) / (γ + 1.0);
    G7 = (γ - 1.0) / 2.0;
    G8 = γ - 1.0;
    # Package these into a single vector
    G = [G1; G2; G3; G4; G5; G6; G7; G8];

    # Compute the speeds of sound
    CL = speedOfSound(QL, γ);
    CR = speedOfSound(QR, γ);

    # Check pressure positivity condition
    if (G4 * (CL + CR) <= UR - UL)
        error("Pressure positivity condition violated")
    end

    # Compute exact solution for pressure and velocity in the star region
    WL = [DL; UL; PL]; # left primitive variables
    WR = [DR; UR; PR]; # right primitive variables
    PM, UM = newton(WL, WR, CL, CR, G);

    # Sample the solution
    if (ξ <= UM) # Left of the contact discontinuity
        if (PM <= PL) # Left rarefaction
            SHL = UL - CL;
            if (ξ <= SHL) # Left data state
                ρ = QL[1];
                ρu = QL[2]; 
                ρE = QL[3];
            else
                CML = CL*(PM/PL)^(G[1]);
                STL = UM - CML;
                if (ξ > STL) # Left star region
                    ρ = DL*(PM/PL)^(1.0/γ);
                    ρu = ρ*UM;
                    WM = [ρ; UM; PM];
                    e = internalEnergy(WM, γ);
                    ρE = ρ * (e + 0.5 * UM^2);
                else # Inside rarefaction fan
                    C = G[5]*(CL + G[7]*(UL - ξ));
                    ρ = DL*(C/CL)^(G[4]);
                    U = G[5]*(CL + G[7]*UL + ξ);
                    ρu = ρ*U;
                    P = PL*(C/CL)^(G[3]);
                    W = [ρ; U; P];
                    e = internalEnergy(W, γ);
                    ρE = ρ * (e + 0.5 * U^2);
                end
            end
        else # Left shock
            PML = PM/PL;
            SL = UL - CL*sqrt(G[2]*PML + G[1]);
            if (ξ <= SL) # Left data state
                ρ = QL[1];
                ρu = QL[2]; 
                ρE = QL[3];
            else # Left star region
                ρ = DL*(PML + G[6])/(PML*G[6] + 1.0);
                ρu = ρ*UM;
                WM = [ρ; UM; PM];
                e = internalEnergy(WM, γ);
                ρE = ρ * (e + 0.5 * UM^2);
            end
        end
    else # right of contact discontinuity
        if (PM > PR) # Right shock
            PMR = PM/PR;
            SR = UR + CR*sqrt(G[2]*PMR + G[1]);
            if (ξ >= SR) # Right data state
                ρ = QR[1];
                ρu = QR[2]; 
                ρE = QR[3];
            else # Right star region
                ρ = DR*(PMR + G[6])/(PMR*G[6] + 1.0);
                ρu = ρ*UM;
                WM = [ρ; UM; PM];
                e = internalEnergy(WM, γ);
                ρE = ρ * (e + 0.5 * UM^2);
            end
        else # Right rarefaction
            SHR = UR + CR;
            if (ξ >= SHR) # Right data state
                ρ = QR[1];
                ρu = QR[2]; 
                ρE = QR[3];
            else
                CMR = CR*(PM/PR)^(G[1]);
                STR = UM + CMR;
                if (ξ <= STR) # Right star region
                    ρ = DR*(PM/PR)^(1.0/γ);
                    ρu = ρ*UM;
                    WM = [ρ; UM; PM];
                    e = internalEnergy(WM, γ);
                    ρE = ρ * (e + 0.5 * UM^2);
                else # Inside rarefaction fan
                    C = G[5]*(CR - G[7]*(UR - ξ));
                    ρ = DR*(C/CR)^(G[4]);
                    U = G[5]*(-CR + G[7]*UR + ξ);
                    ρu = ρ*U;
                    P = PR*(C/CR)^(G[3]);
                    W = [ρ; U; P];
                    e = internalEnergy(W, γ);
                    ρE = ρ * (e + 0.5 * U^2);
                end
            end
        end
    end

    # Return vector of conserved variables
    return [ρ; ρu; ρE]

end

