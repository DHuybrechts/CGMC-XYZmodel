function A = GetOperators( o , N)
%-------------------------------------------------------------------------%
%Calculates the operators on the product wave function for every site, 
%given the operator 'o' for that site.
%Parameters:
    %o      operator
    %N      number of sites in cluster
%-------------------------------------------------------------------------%

    A = zeros(2^N,2^N,N);
    C = cell(1, N);
    for i = 1:N
        C(i) = {eye(2)};
    end
    
    %apply kronecker product:
    for i = 1:N
        D = C;
        D(i) = {o};
        A(:,:,i) = superkron(D);
    end
end