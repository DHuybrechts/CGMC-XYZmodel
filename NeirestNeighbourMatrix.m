function NNM = NeirestNeighbourMatrix( out )
%-------------------------------------------------------------------------%
%   Gives a matrix where each row and column represents the site in the
%   lattice. Each row vector contains a 1 in the column of a neirest neighbour
%   of that row's site outside the cluster, and a 0 for all the others.
%Parameters:
%   out     neighbours of site i outside cluster k.
%-------------------------------------------------------------------------%
    NNM = zeros(size(out,1));
    
    for i = 1:size(out,1)
        n_out = nonzeros(out(i,:));
        for j = 1:length(n_out)
            NNM(i,n_out(j)) = NNM(i,n_out(j)) + 1;
        end
    end

    NNM = sparse(NNM);
end