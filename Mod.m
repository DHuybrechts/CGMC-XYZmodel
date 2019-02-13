function out = Mod( a, b, offset )
%-------------------------------------------------------------------------%
%   calculates the modulus of the number a with specific offset
%Parameters:
%   a           numerator
%   b           denominator
%   offset      offset of the modulus.
%-------------------------------------------------------------------------%
r=mod(a,b);
while(r<offset)
    r=r+b;
end
out=r;
end