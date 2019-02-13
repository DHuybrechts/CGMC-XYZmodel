function [all, in, out] = CalcClusterNeighbours1D( Nxc, Nyc, Nx, Ny )
%-------------------------------------------------------------------------%
%   Calculate the neighbours of a site 'i' on the Nx x Ny lattice. With i 
%   located in cluster k. The sites are numbered from 1 -> Nx.Ny starting
%   from the first row and first column. (Periodic boundary conditions)
%Parameters:
%   Nxc             number of sites in x direction of cluster
%   Nyc             number of sites in y direction of cluster
%   Nx              number of rows
%   Ny              number of columns
%Output:
%   all             all neighbours of site i.
%   in              neighbours of site i inside cluster k.
%   out             neighbours of site i outside cluster k.
%Mapping between lineair index 'a' in row vector and the sites on the 2D xy
%plain are given by a=Ny(x-1)+y. where x=ceil(j/Ny) and y=Mod(j,Ny,1)
%-------------------------------------------------------------------------%
    t = zeros(Nx*Ny,2);
    for j=1:Nx*Ny
        t(j,1)=Ny*(ceil(j/Ny)-1)+Mod(Mod(j,Ny,1)-1,Ny,1);   %Left neighbour
        t(j,2)=Ny*(ceil(j/Ny)-1)+Mod(Mod(j,Ny,1)+1,Ny,1);   %Right neighbour
    end
    
    %open boundary conditions
%     t(1,1) = 0;
%     t(end,2) = 0;
    
    all = t;
    vc = ClusterConfiguration(Nxc, Nyc, Nx, Ny );
    r = zeros(size(all));
    for i = 1:size(all, 1)
        [x, ~] = find(vc == i);                     %Index of cluster from site i.
        r(i,:) = ismember(all(i,:),vc(x,:));
    end
    in = all.*r;
    out = all - in;
end