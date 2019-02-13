function Sigxyz = CalcExpSig( C, switch_xy, magn_y, magn_z, clustconfig)
%-------------------------------------------------------------------------%
%   Calculate the expectation value of the x,y,z magnetization of each
%   site.
%Parameters:
%   C               Coefficient value the ODE will calculate
%   switch_xy,
%   magn_y, magn_z  Alternative way to calculate expectation values, uses
%                   the permutation properties of the pauli matrices.
%   clustconfig     matrix where the row index represents the cluster and
%                   contains the indices of the sites in that respective 
%                   cluster.
%-------------------------------------------------------------------------%
    n_sc = size(clustconfig, 2);            %number of sites in cluster
    n_coeff = 2^n_sc;                       %number of coefficients in cluster
    n_clust = size(clustconfig, 1);         %number of clusters
    sigxyz = zeros(n_clust*n_sc, 1);
    
    for i = 1:n_clust                       %Sum over clusters
        range = n_coeff*(i-1)+1:n_coeff*(i-1)+1+(n_coeff-1);
        d = C(range);
        cd = d';
        for j = 1:n_sc                      %Sum over sites in cluster
            ds = d(switch_xy{j});
            Ind = clustconfig(i,j);
            
            sigxyz(Ind, 1) = cd*ds;                         %d'*sig(:,:,j,1)*d;

            sigxyz(Ind, 2) = cd*(ds.*magn_y{j});            %d'*sig(:,:,j,2)*d;

            sigxyz(Ind, 3) = cd*(d.*magn_z{j});             %d'*sig(:,:,j,3)*d;
        end
    end
    Sigxyz = sigxyz;
end