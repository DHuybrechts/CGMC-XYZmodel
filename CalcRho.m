function Rho = CalcRho( Ct, Nxc, Nyc, Nx, Ny)
%-------------------------------------------------------------------------%
%   Calculates the density matrix of a wave function at a certain time.
%Parameters:
%   Ct              The wave function in coëfficiënt notation for a time..
%   Nxc             number of sites in x direction of cluster
%   Nyc             number of sites in y direction of cluster
%   Nx              number of rows
%   Ny              number of columns
%-------------------------------------------------------------------------%
    N = Nx*Ny;                                                              %number of sites in lattice
    Nc = Nxc*Nyc;                                                           %number of sites in cluster
    NC = N/Nc;                                                              %number of clusters in lattice
    N_coeff = 2^Nc;                                                         %number of coefficients in cluster wave function
    
    r = cell(1, NC);
    
    for i = 1:NC
        range = N_coeff*(i-1)+1:N_coeff*(i-1)+1+(N_coeff-1);
        r{i} = Ct(range)*Ct(range)';                                        %outer product
    end
    
    Rho = r;
end