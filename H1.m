function OA = H1( gamma, Jx, Jy, Jz, Nxc, Nyc, Nx, Ny )
%-------------------------------------------------------------------------%
%   Calculate the first part of the Hamiltonian, the time independent part,
%   reducing it to single operators (matrices) of each cluster acting on 
%   their respective cluster wave functions.
%Parameters:
%   gamma           system parameter
%   Jx,Jy,Jz        system parameter
%   Nxc             number of sites in x direction of cluster
%   Nyc             number of sites in y direction of cluster
%   Nx              number of rows
%   Ny              number of columns
%-------------------------------------------------------------------------%
    %Some system parameters:
    c = ClusterConfiguration(Nxc, Nyc, Nx, Ny);
    SiteI = SiteClusterIndex(Nxc, Nyc, Nx, Ny);
    if Nxc == 1 && Nx == 1
        [~, in, ~] = CalcClusterNeighbours1D(Nxc, Nyc, Nx, Ny);
    else
        [~, in, ~] = CalcClusterNeighbours(Nxc, Nyc, Nx, Ny);
    end
    sig = GetAllOperators(Nxc*Nyc);
    sig = nDMatrixToCellArraySparse(sig);
    
    %The calculation:
    OA = cell(1,size(c,1));
    for cluster=1:size(c,1)
        A = zeros(2^(Nxc*Nyc),2^(Nxc*Nyc));
        for i_clust = 1:length(c(cluster,:))                                    %index of the site inside the cluster.
            i_site = c(cluster,i_clust);                                        %index of the site in the lattice.
            A = A + (-gamma/2)*sig{i_clust,4};

            n_in = nonzeros(in(i_site,:));                                      %neighbours of i_site in cluster. !in lattice picture!
            for j = 1:length(n_in)
                ind = SiteI(n_in(j));                                           %Cluster index of the site n_in(j);

                A = A + (-1i)*0.5*(Jx*sig{i_clust,1}*sig{ind,1}+...
                    Jy*sig{i_clust,2}*sig{ind,2}+...
                    Jz*sig{i_clust,3}*sig{ind,3});
            end
        end
        OA{cluster} = sparse(A);
    end
end