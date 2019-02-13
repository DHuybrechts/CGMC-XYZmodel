function OA = H2( gamma, Jx, Jy, Jz, clustconfig, NNM, sig, Sigxyz )
%-------------------------------------------------------------------------%
%   Calculate the second part of the Hamiltonian depending on the time 
%   dependent coefficients, reducing it to single operators (matrices) 
%   of each cluster acting on their respective cluster
%   wave functions.
%Parameters:
%   gamma           system parameter
%   Jx,Jy,Jz        system parameter
%   clustconfig     matrix where the row index represents the cluster and
%                   contains the indices of the sites in that respective cluster.
%   NNM             Connection matrix, each row is a site on the lattice 
%                   and the columns with a 1 are the neirest neighbours 
%                   outside the cluster.
%   sig             contains the matrix operators for sigma x, y and z.
%   Sigxyz          x,y,z magnetization for each site.
%-------------------------------------------------------------------------%
    n_sc = size(clustconfig, 2);            %number of sites in cluster
    n_coeff = 2^n_sc;                       %number of coefficients in cluster
    n_clust = size(clustconfig, 1);         %number of clusters

    OA = cell(1, n_clust);
    for i = 1:n_clust
        OA{i} = zeros(n_coeff,n_coeff);
    end
    
    sx = Jx*NNM*Sigxyz(:,1);                                                %x-magnetization mean-field contribution on each site
    sy = Jy*NNM*Sigxyz(:,2);                                                %y-magnetization mean-field contribution on each site
    sz = Jz*NNM*Sigxyz(:,3);                                                %z-magnetization mean-field contribution on each site
    for cluster=1:n_clust
        for i_clust = 1:n_sc                                                %index of the site inside the cluster.
            i_site = clustconfig(cluster,i_clust);                          %index of the site in the lattice.

            OA{cluster} = OA{cluster} + (sig{i_clust,1}*sx(i_site) +...
                sig{i_clust,2}*sy(i_site)+...
                sig{i_clust,3}*sz(i_site));
        end
        OA{cluster} = (-1i)*OA{cluster};
    end
end