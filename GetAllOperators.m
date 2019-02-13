function B = GetAllOperators( N )
%-------------------------------------------------------------------------%
%Calculates the representations of sigma^x, sigma^y, sigma^z, sigma^pm,
%sigma^m and the products sigma^xi*sigma^xj, sigma^yi*sigma^yj and
%sigma^zi*sigma^zj for the different sites in the cluster.
%Parameters:
%   N   number of sites in cluster
%-------------------------------------------------------------------------%
    sx = [ 0 1; 1 0];                   %x pauli matrix
    sy = [ 0 -1i; 1i 0];                %y pauli matrix
    sz = [ 1 0; 0 -1];                  %z pauli matrix
    spm = [ 1 0; 0 0];                  %raising operator * lowering operator
    sm = [ 0 0; 1 0];                   %lowering operator
    B = zeros(2^N,2^N,N,5);
    
    %representations for the different sites:
    B(:,:,:,1) = GetOperators(sx, N);
    B(:,:,:,2) = GetOperators(sy, N);
    B(:,:,:,3) = GetOperators(sz, N);
    B(:,:,:,4) = GetOperators(spm, N);
    B(:,:,:,5) = GetOperators(sm, N);
end

