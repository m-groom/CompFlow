# Functions for computing iterative solution to the exact Riemann problem
# Note: currently only valid for fixed γ and ideal gas EOS

# Newton iteration to compute p* and u*
function newton(WL, WR, CL, CR, G)
    # Extract velocities
    UL = WL[2]; UR = WR[2]; 

    # Newton iteration to compute p* and u*
    PS, _ = ANRS(WL, WR, CL, CR, G);   # Initial guess for p*
    ϵ = 1e-12;                      # tolerance
    maxIter = 100;                  # maximum number of iterations

    # Compute pressure in the star region
    for i=1:maxIter
        g, dg = f(PS, WL, WR, CL, CR, G);  # function to solve and its derivative
        res = abs(g); 
        if(res < ϵ)        
            break # we have found the solution so we exit
        end
        Δp = -g / dg; # Newton step   
        # line search globalization 
        δ = 1; # scale factor 0<δ<=1 
        for ii=1:20   
            gk, _ = f(PS + δ * Δp, WL, WR, CL, CR, G);       
            if(abs(gk) >= res)
                # if the residual of the next iteration increases then the Newton step is reduced by 2
                δ = 0.5 * δ; 
            else 
                # residual does not increase => exit inner loop
                break
            end
        end
        # Update the solution. For δ=1, globalization method reduces to standard Newton. 
        PS = PS + δ * Δp; 
    end

    # Compute velocity in the star region
    fR, _ = fK(PS, WR, CR, G); fL, _ = fK(PS, WL, CL, G);
    US = 0.5 * ((UL + UR) + (fR - fL));

    return PS, US

end

# Function gk used in the root finding algorithm
function f(PS, WL, WR, CL, CR, G)
    # Extract velocities
    UL = WL[2]; UR = WR[2]; 
    # Evaluate fR, fL and their derivatives
    fR, dfR = fK(PS, WR, CR, G);
    fL, dfL = fK(PS, WL, CL, G);
    # Compute g and its derivative
    g = fR + fL + UR - UL;
    dg = dfR + dfL;

    return g, dg

end

# Function fK containing the Rankine-Hugoniot relation, Riemann invariants and Lax entropy condition
function fK(P, WK, CK, G)
    # Extract primitive variables
    DK = WK[1]; PK = WK[3];
    AK = G[5] / DK;
    BK = G[6] / DK;

    # Lax entropy condition
    if (P > PK) 
        # shock (Rankine-Hugoniot conditions)
        ϕ = (P - PK) * sqrt(AK / (P + BK));
        dϕ = (1.0 - 0.5*(P - PK)/(BK + P)) * sqrt(AK / (P + BK));
    else
        # rarefaction (Riemann invariants) 
        ϕ = G[4] * CK * ((P / PK)^(G[1]) - 1.0);
        dϕ = (1.0 / (DK*CK)) * (P / PK)^(-G[2]);
    end

    return ϕ, dϕ

end

# Adaptive Noniterative Riemann solver (from Sec. 9.5.2 of Toro)
# Returns pressure in the star region by selecting from PVRS, TRRS or TSRS
function ANRS(WL, WR, CL, CR, G)
    # Extract primitive variables
    DL = WL[1]; UL = WL[2]; PL = WL[3];
    DR = WR[1]; UR = WR[2]; PR = WR[3];
    # User specified pressure ratio
    Quser = 2.0;
    # Compute initial guess for pressure in the star region from PVRS
    CUP = 0.25 * (DL + DR) * (CL + CR);
    PPV = max(0.0, 0.5 * (PL + PR) + 0.5 * (UL - UR) * CUP);
    Pmin = min(PL, PR);
    Pmax = max(PL, PR);
    Qmax = Pmax / Pmin;

    if (Qmax <= Quser && (Pmin <= PPV && PPV <= Pmax))
        # Use solution from PVRS Riemann solver
        PM = PPV;
        UM = 0.5 * (UL + UR) + 0.5 * (PL - PR) / CUP;
    else
        if (PPV < Pmin)
            # Use solution from the Two-Rarefaction Riemann solver
            PQ = (PL / PR)^(G[1]);
            UM = (PQ*UL/CL + UR/CR + G[4]*(PQ - 1.0))/(PQ/CL + 1.0/CR);
            PTL = 1.0 + G[7]*(UL - UM)/CL;
            PTR = 1.0 + G[7]*(UM - UR)/CR;
            PM = 0.5 * (PL*PTL^(G[3]) + PR*PTR^(G[3]));
        else
            # Use solution from the Two-Shock Riemann solver
            GEL = sqrt(G[5]/(DL*(G[6]*PL + PPV)));
            GER = sqrt(G[5]/(DR*(G[6]*PR + PPV)));
            PM = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER);
            UM = 0.5*(UL + UR) + 0.5*(GER*(PM - PR) - GEL*(PM - PL))
        end
    end

    return PM, UM

end