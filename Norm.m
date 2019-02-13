function n = Norm( c , T)
%-------------------------------------------------------------------------%
%   Calculates the norm of the product of wave functions at every timestep.
%Parameters:
%       c       wave function in coefficients (column vector)
%       T       matrix which makes the sum over the moduli of the
%               coefficents of the different cluster wave functions
%-------------------------------------------------------------------------%
    n = prod(T*(conj(c).*c));
end
