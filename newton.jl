# Newton iteration to compute p* and u*
function newton(WL, WR, CL, CR, G)
    # Extract primitive variables
    DL = WL[1]; UL = WL[2]; PL = WL[3];
    DR = WR[1]; UR = WR[2]; PR = WR[3];
    # Initial guess for p
    PS = ANRS(WL, WR, CL, CR, G);
    tol = 1e-12;        # tolerance
    maxIter = 100;      # maximum number of iterations
    for i=1:maxIter
        # gk = f(hs,hL,hR,uL,uR,g);  # function to solve
        # res = abs(gk); 
        # if(res<tol)        
        #     break # we have found the solution so we exit
        # end
        # dg = df(hs,hL,hR,uL,uR,g); # derivative of the function
        # dh = -gk/dg; # Newton step   
        # # line search globalization 
        # delta = 1; # scale factor 0<delta<=1 
        # for ii=1:20         
        #     if(abs(f(hs+dh*delta,hL,hR,uL,uR,g)) >= res)
        #         # if the residual of the next iteration increases then the Newton step is reduced by 2
        #         delta = 0.5*delta; 
        #     else 
        #         # residual does not increase => exit inner loop
        #         break
        #     end
        # end
        # # Update the solution. For delta=1, globalization method reduces to standard Newton. 
        # hs = hs + delta*dh; 
    end
    return PS, US
end

# # Function gk used in the root finding algorithm
# function f(hs,hL,hR,uL,uR,g)
#     return phi(hs,hL,g) + phi(hs,hR,g) + uR - uL;
# end

# # Function phi containing the Rankine-Hugoniot relation, Riemann invariants and Lax entropy condition
# function phi(hs,hLR,g)
#     # Lax entropy condition
#     if (hs>hLR) 
#         # shock (Rankine-Hugoniot conditions)
#         y = sqrt(0.5*g*(hs+hLR) / (hs*hLR)) * (hs-hLR);  
#     else
#         # rarefaction (Riemann invariants) 
#         y = 2*sqrt(g)*(sqrt(hs) - sqrt(hLR)); 
#     end
#     return y
# end

# # Derivative of the function gk
# function df(hs,hL,hR,uL,uR,g)
#     eps = 1e-6;
#     # Central finite difference approximation
#     return (f(hs+eps,hL,hR,uL,uR,g)-f(hs-eps,hL,hR,uL,uR,g))/(2*eps); 
# end

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
    PPV = max(0.0, 0.5 * (PL + PR) - 0.5 * (UR - UL) * CUP);
    Pmin = min(PL, PR);
    Pmax = max(PL, PR);
    Qmax = Pmax / Pmin;

    if (Qmax <= Quser && (Pmin <= PPV && PPV <= Pmax))
        # Use solution from PVRS Riemann solver
        PM = PPV;
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
            GEL = sqrt((G[5]/DL)/(G[6]*PL + PPV));
            GER = sqrt((G[5]/DR)/(G[6]*PR + PPV));
            PM = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER);
        end
    end

    return PM

end