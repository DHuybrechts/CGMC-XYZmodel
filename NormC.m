  function n = NormC( c , Nc)
%-------------------------------------------------------------------------%
%   Calculates the norm factor for every coëfficiënt in the wave function.
%Parameters:
%       c       wave function in coefficients (row or column vector)
%       Nc      number of clusters
%-------------------------------------------------------------------------%
    n_coeff = length(c)/Nc;             %Number of coefficients in cluster wave function.
    
    n_wf = zeros(size(c));              %Will contain the norm of the individual coefficients of the cluster wave functions.
    C = conj(c).*c;
    
    for i =1:Nc
        range = n_coeff*(i-1)+1:n_coeff*(i-1)+1+(n_coeff-1);
        n_wf(range) = sum(C(range));
    end
    n = sqrt(n_wf);                     %contains the normalization factor for each coëfficiënt.
end