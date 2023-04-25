# CompFlow
Solvers based on the Godunov finite-volume method for compressible flow problems.

---

Boundary conditions:
- 0 = periodic (default)
- 1 = extrapolate
- 2 = reflective
- 3 = wall
- 4 = user defined

The solver computes the approximate solution to the 1D compressible Euler equations
$$\frac{\partial \boldsymbol Q}{\partial t} + \frac{\partial \boldsymbol F}{\partial x} = 0$$
where $\boldsymbol Q = [\rho, \rho u, \rho E]^T$ is the vector of conserved variables and the flux vector $\boldsymbol F$ is given by $$\boldsymbol F = \left[ \rho u, \rho u^2 + p, u(\rho E + p) \right]^T$$

Here $\rho$ denotes the density, $u$ is the velocity in the $x$ direction, $p$ is the pressure, and $E=e+k$ is the total energy. The kinetic energy is given by $k=\frac{1}{2} u^2$ and the internal energy $e$ is calculated through the equation of state.