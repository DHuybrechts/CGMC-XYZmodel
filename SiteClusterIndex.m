function SiteI = SiteClusterIndex( Nxc, Nyc, Nx, Ny )
%-------------------------------------------------------------------------%
%   Gives a vector of which the row index is the index of the site, the
%   corresponding value in the vector is the index of the respective site 
%   in his cluster.
%Parameters:
%   Nxc             number of sites in x direction of cluster
%   Nyc             number of sites in y direction of cluster
%   Nx              number of rows
%   Ny              number of columns
%-------------------------------------------------------------------------%

    clustconfig = ClusterConfiguration(Nxc, Nyc, Nx, Ny );
    SiteI = zeros(Nx*Ny,1);
    for i=1:Nx*Ny
        [~, SiteI(i)] = find(clustconfig == i);
    end
end

