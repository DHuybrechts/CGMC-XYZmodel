function P = ChanceInterval( cin, clustconfig, sig, T )
%-------------------------------------------------------------------------%
%   Calculates the chance intervals of the sites to make a jump.
%Parameters:
%   cin     input wavefunction in coefficient notation
%   Nx      number of rows
%   Ny      number of columns
%   NC      number of clusters in lattice
%-------------------------------------------------------------------------%
    
    n = numel(clustconfig);
    J = zeros(n,1);
    J2 = J;
    for i = 1:n
        J(i) = Norm(Jump(cin, clustconfig, sig, i), T);
        J2(i) = sum(J(1:i));
    end
    P = J2/sum(J);
    
    
    %Give possible warning and errors:
    if (abs(P(end) - 1) > 1e-3)
        warning('Sum of the chances deviates more than 1e-3 from 1.')
    elseif isnan(P(end))
        error('Chance interval contains NaN data type.')
    end
end

