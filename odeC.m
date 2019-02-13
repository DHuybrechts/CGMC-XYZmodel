function dC = odeC(t, C, gamma, Jx, Jy, Jz, Nxc, Nyc, Nx, Ny, clustconfig, NNM, sig, switch_xy, magn_y, magn_z, A1)
%-------------------------------------------------------------------------%
%   system of differential equations for the coefficients of the individual
%   site wave functions.
%Parameters:
%   t           time variable
%   C           coefficient value the ODE will calculate (column vector).
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
%-------------------------------------------------------------------------%
    Nc = Nxc*Nyc;   %number of sites in cluster.
    N=Nx*Ny/(Nc);   %number of clusters in lattice.  
    
    if (Nxc == Nx) && (Nyc == Ny)%exact solution
        A2 = cell(1);
        A2{1} = 0;
    else%cluster solution
        Sigxyz = CalcExpSig( C, switch_xy, magn_y, magn_z, clustconfig );       %We calculate the expectation value of sigma_i^r, with r= x, y, z.
        A2 = H2(gamma, Jx, Jy, Jz, clustconfig, NNM, sig, Sigxyz );             %Time dependent part of the Hamiltonian.
    end
                                                                            
    dC = zeros(N*2^Nc,1);
    for j=1:N      %sum over the clusters.
        range = (2^Nc)*(j-1)+1:(2^Nc)*(j-1)+1+(2^Nc-1);
        if (Nxc == Nx) && (Nyc == Ny)
            dC(range) = A1{j}*C(range);
        else
            dC(range) = (A1{j}+A2{j})*C(range);
        end
    end
end   