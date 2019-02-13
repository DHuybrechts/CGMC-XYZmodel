function [ switch_xy, magn_y, magn_z ] = CalcAlternativeOperator(sig)
%-------------------------------------------------------------------------%
%   Calculates vectors which are used for an alternative matrix operator
%   product, based on the permutation properties, which is computationally
%   more efficient in certain cases.
%Parameters:
%   sig         contains the matrix operators for sigma x, y and z.
%-------------------------------------------------------------------------%
    %Data structures:
    Index = transpose(1: size(sig{1,1},1));
    switch_xy = cell(size(sig,3));
    magn_y = switch_xy;
    magn_z = switch_xy;
    
    %calculate alternative 'operators':
    for i = 1:size(sig,1)
        switch_xy{i} = sig{i,1}*Index;
        magn_y{i} = sign(sig{i,2}*Index);
        magn_z{i} = sign(sig{i,3}*Index);
    end
end

