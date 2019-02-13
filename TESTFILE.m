%TESTFILE
clear all

Nxc = 2;            %x dimension of cluster
Nyc = 2;            %y dimension of cluster
Nx = 4;             %x dimension of lattice
Ny = 4;             %y dimension of lattice

gamma = 1;          %dissipation rate
Jx = 0.9*gamma;     %coupling parameter x direction
Jy = 1.25*gamma;    %coupling parameter y direction
Jz = gamma;         %coupling parameter z direction
T = 10000;             %time window of the trajectory
dt = 1;             %after each dt 'save' the state of the system

%calculate the trajectory and the steady-state spin structure factor 'Ss'
%at each point in the time window. 
%time       contains the time steps
%Ss         contains the steady-state spin structure factor at each time
%Wf         contains the wave function at each time

[ time, Ss, Wf] = CalculateSssxx( gamma, Jx, Jy, Jz, Nxc, Nyc, Nx, Ny, T, dt);

