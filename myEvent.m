function [ value, isterminal, direction ] = myEvent( t, C, eps, T)
%-------------------------------------------------------------------------%
%If value becomes zero an event occurs and our system will jump at this 
%time.
%   t           time
%   C           solution of odeC
%   eps         random number
%   T           matrix which makes the sum over the moduli of the
%               coefficents of the different cluster wave functions
%-------------------------------------------------------------------------%
    value = Norm(C, T) - eps;
    isterminal = 1;
    direction = -1;
end