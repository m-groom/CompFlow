# Define the flux as a function of Q
include("equation_of_state.jl")
function F(Q, gamma)
    flux = zeros(3,1);
    rho = Q[1];
    if (rho <= 0)
        rho = 0;
        u = 0;
        p = 0;
        rhoE = 0;
    else
        u = Q[2] / Q[1];
        p = eos(Q, gamma);
        rhoE = Q[3];
    end
    flux[1] = rho * u; # mass flux
    flux[2] = rho * u^2 + p; # momentum flux
    flux[3] = u * (rhoE + p); # energy flux
    return flux
end