# Define the flux as a function of Q
include("equation_of_state.jl")
function F(Q, γ = 5/3)
    flux = zeros(3,1);
    ρ = Q[1];
    if (ρ <= 0)
        ρ = 0;
        u = 0;
        p = 0;
        ρE = 0;
    else
        u = Q[2] / Q[1];
        p = eos(Q, γ);
        ρE = Q[3];
    end
    flux[1] = ρ * u;          # mass flux
    flux[2] = ρ * u^2 + p;    # momentum flux
    flux[3] = u * (ρE + p);   # energy flux
    return flux
end