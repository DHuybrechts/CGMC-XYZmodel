function t = GetSite( P, rnd )
%-------------------------------------------------------------------------%
%   Gives the site that jumps for given random number.
%Parameters:
%   P       chance interval of the sites, gives probability to jump.
%   rnd     random number for the jump location.
%-------------------------------------------------------------------------%
s=find(sign(P-rnd)>0); 
t=s(1);
end

