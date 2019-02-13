function Sss = CalcSss( C, sig, clustconfig, N )
%-------------------------------------------------------------------------%
%   Calculates the steady-state spin structure factor for a certain time.
%Parameters:
%   C              The wave function in coëfficiënt notation for a time..
%   sig             contains the matrix operators for sigma x, y and z.
%   clustconfig     matrix where the row index represents the cluster and
%                   contains the indices of the sites in that cluster.
%   N               Number of sites.
%-------------------------------------------------------------------------%
    s=0;
    r = sig{1,1};
    n_coeff = size(r,1);
    
    for i = 1:N
        for j= 1:N
            if j ~= i
                [c1, s1] = find(clustconfig == i);
                [c2, s2] = find(clustconfig == j);
                
                if c1 == c2 %same cluster
                    range = n_coeff*(c1-1)+1:n_coeff*(c1-1)+1+(n_coeff-1);
                    s = s + C(range)'*sig{s1,1}*sig{s2,1}*C(range);
                else        %different cluster
                    range1 = n_coeff*(c1-1)+1:n_coeff*(c1-1)+1+(n_coeff-1);
                    range2 = n_coeff*(c2-1)+1:n_coeff*(c2-1)+1+(n_coeff-1);
                    s = s + (C(range1)'*sig{s1,1}*C(range1))*(C(range2)'*sig{s2,1}*C(range2));
                end
                
            end
        end
    end
    
    Sss = s/(N*(N-1));
end