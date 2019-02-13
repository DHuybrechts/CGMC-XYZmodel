function [time, Ct, tj, Cj, ij] = EffH( gamma, Jx, Jy, Jz, Nxc, Nyc, Nx, Ny, clustconfig, NNM, sig, switch_xy, magn_y, magn_z, A1, Cin, T, Tspan, eps )
%-------------------------------------------------------------------------%
%Function EffH 
%   Calculates the time evolution of the wave function.
%Parameters:
%   gamma           system parameter
%   Jx,Jy,Jz        system parameter
%   Nxc             number of sites in x direction of cluster
%   Nyc             number of sites in y direction of cluster
%   Nx              number of rows
%   Ny              number of columns
%   clustconfig     matrix where the row index represents the cluster and
%                   contains the indices of the sites in that respective cluster.
%   NNM             Connection matrix, each row is a site on the lattice 
%                   and the columns with a 1 are the neirest neighbours 
%                   outside the cluster.
%   sig             contains the matrix operators for sigma x, y and z.
%   switch_xy,
%   magn_y, magn_z  Alternative way to calculate expectation values, uses
%                   the permutation properties of the pauli matrices.
%   A1              matrix containing the Hamiltonian independt of
%                   time, per cluster.
%   Cin             input wave function.
%   T               matrix which makes the sum over the moduli of the
%                   coefficents of the different cluster wave functions
%   Tspan           Time evolution up till time T
%   eps             the random value to determine the jump time
%-------------------------------------------------------------------------%
options = odeset('events', @(t,C) myEvent(t, C, eps, T));
[time, Ct, tj, Cj, ij] = ode45(@(t, C) odeC(t, C, gamma, Jx, Jy, Jz, Nxc, Nyc, Nx, Ny, clustconfig, NNM, sig, switch_xy, magn_y, magn_z, A1), Tspan, Cin, options);
end