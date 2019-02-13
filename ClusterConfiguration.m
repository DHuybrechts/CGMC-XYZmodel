function A = ClusterConfiguration( Nxc, Nyc, Nx, Ny )
%-------------------------------------------------------------------------%
%   Calculates the cluster configuration on the Nx x Ny lattice of 
%   Nxc x Nyc clusters. The sites are numbered from 1 -> Nx.Ny starting 
%   from the first row and first column. 
%Parameters:
%   Nxc             number of sites in x direction of cluster
%   Nyc             number of sites in y direction of cluster
%   Nx              number of rows
%   Ny              number of columns
%-------------------------------------------------------------------------%
    Nc = (Nx*Ny)/(Nxc*Nyc);
    r = 1;
    d = zeros(Nx,Ny);
    for i =1:Nx
        for j=1:Ny
            d(i,j) = r;
            r = r + 1;
        end
    end
    c = zeros(Nc, Nxc*Nyc);
    r = 1;
    for i = 1:Nx/Nxc
        for j = 1:Ny/Nyc
            x1 = (i-1)*Nxc+1;
            x2 = x1+Nxc-1;
            y1 = (j-1)*Nyc+1;
            y2 = y1+Nyc-1;
            
            clm = d(x1:x2,y1:y2);   
            clm = clm';
            c(r,:)=clm(:)';
            r = r + 1;
        end
    end
    
    A = c;
end