function Ms= nDMatrixToCellArraySparse( M )
%-------------------------------------------------------------------------%
%   transforms a n-dimensional matrix into a sparse array
%Parameters:
%   M   n-dimensional matrix
%-------------------------------------------------------------------------%
    s = size(M);
    Ms = cell(s(3), s(4));
    
    for i = 1:s(4)
        for j = 1:s(3)
            Ms{j, i} = sparse(M(:,:,j,i));
        end
    end
end

