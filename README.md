# This repository contains methods to solve ODEs and PDEs in MATLAB, descriptions of each method are shown below:

- wave1d:                 The finite-difference method for the 1D wave equation involves locating a matrix U that approximates the PDE.
- laplace2d:              We will use finite difference equations with insulated boundaries to approximate the solution to Laplace's equation in 2D and 3D. 
                          This will result in a system of equations that we can solve.
- Euler’s method:         Numerical method for approximating the solution to a 1st order IVP. Then we will add the trapezoidal rule to improve our code.
- Crank_nicolson1d:       The Crank-Nicolson method finds a matrix U that approximates the PDE using finite-difference to approximate solutions to the 
                          heat-conduction/diffusion equation.
- Boundary Value Problem: Approximating numerical solutions by creating a matrix in a BVP using linear ordinary differential equations (ODEs).
- richardson22:           The Richardson extrapolation helps us find an approximation of the derivative of a function at a specific point using multiple
                          approximations of the derivative. This method gives us closer values to the real values with digits of error.
