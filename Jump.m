function c = Jump( cin, clustconfig, sig, i )
%-------------------------------------------------------------------------%
%   Gives the wave function when a jump occured on site i
%Parameters:
%   cin             wave function in coefficient notation.
%   Psi             Cluster wave function from cluster containing site i,
%                   containing the spin up, spin down states.
%   clustconfig     configuration of the lattice sites in the cluster.
%   i               index of site where the jump occurs.
%-------------------------------------------------------------------------%
n_coeff = 2^size(clustconfig,2);                        %number of coefficients in cluster wave function.
[x, y] = find(clustconfig == i);                        %x = cluster index and y = site index in cluster.
range = n_coeff*(x-1)+1:n_coeff*(x-1)+1+(n_coeff-1);    %Denotes the indices of the cluster wave function in the system wave function.
cin(range) = sig{y,5}*cin(range);                       %Lowering operator on site i works on the wave function of the respective cluster.
c = cin;
end

