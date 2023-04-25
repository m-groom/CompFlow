# Functions for updating the solution at each timestep

# Function for computing the approximate solution
function computeSolution!(Q, x, BCs, γ, CFL, Nmax, Nout, tstart, tend)
    # Set the initial time
    t = tstart;
    # Start explicit timestepping
    for n = 1:Nmax
        # Compute the time step
        Δt = getTimeStep(Q, x, γ, CFL);
        if (t + Δt > tend)
            Δt = tend - t;
        end
        # Stop criterion
        if (t >= tend)
            break
        end
        report("Current time step: $(n)   Current time: $(rpad(string(round(t + Δt, digits=6)), 8, "0"))")
        # Reconstruct the extrapolated values at the cell boundary
        QR, QL = reconstruct(Q, BCs, γ);
        # Evolve the extrapolated values at the cell boundary
        evolve!(QR, QL, x, Δt, γ);
        # Compute the solution at the next timestep
        update!(Q, QR, QL, x, BCs, Δt, γ);
        # Update the time
        t = t + Δt;
        # Write the solution at every Nout time steps
        if (n%Nout == 0)
            report("Saving the solution at time t = $(rpad(string(round(t, digits=6)), 8, "0"))")
            writeSolution(x, Q, "solution_$(rpad(string(round(t, digits=6)), 8, "0")).vtr")
        end
    end

    return t

end

# Function for estimating the timestep size
function getTimeStep(Q, x, γ, CFL)
    # Number of cells
    imax = size(Q, 2); 
    # Estimate the maximum wave speed
    Smax = 0.0;
    Δx = x[end] - x[1];
    for i = 1:imax
        ai = speedOfSound(Q[:, i], γ);
        ui = Q[2, i] / Q[1, i];
        Smax = max(Smax, abs(ui) + ai);
        Δx = min(Δx, x[i+1] - x[i]);
    end
    
    return CFL * Δx / Smax

end

# Function for updating the solution
# TODO: generalise boundary conditions
function update!(Q, QR, QL, x, BCs, Δt, γ)
    imax = length(x) - 1; # number of cells
    for i = 1:imax
        Δx = x[i+1] - x[i]; # grid spacing
        if (i==1)
            Q[:,i] .= Q[:,i]; # Dirichlet boundary condition
        elseif (i==imax)
            Q[:,i] .= Q[:,i]; # Dirichlet boundary condition
        else
            Fp = riemannSolver(QR[:,i], QL[:,i+1], γ) # Flux at x[i+1/2]
            Fm = riemannSolver(QR[:,i-1], QL[:,i], γ) # Flux at x[i-1/2]
            Q[:,i] .= Q[:,i] - Δt/Δx * (Fp - Fm);
        end
    end
end