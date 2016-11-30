



function out_indices = bounds_to_indices(bounds)
%     Takes in a set of function boundaries and returns
%     a set of indices corresponding to those boundaries.
%     
%     bounds is a 2xM matrix where each column represnts a pair of boundaries
%     out_indices is a NxM matrix where each column contains the set of indicies
%         specified by the bounds matrix.
%         
%     Example
%     bounds = [1 3; 5 7; 6 8]';
%     out_indices = [1 2 3; 5 6 7; 6 7 8]';

    N_indices=diff(bounds,1)+1;
    
    %Check that for each pair of bounds specified, the number of indices is
    %constant. This is required to ensure the matrix is squre.
    if sum(diff(N_indices,2)) ~= 0
        fprintf('The spacing between all boundaries must be equal.\n');
        return
    end
    
    N_indices = mode(N_indices);        % All N_indices elements should be the same
    N_bounds = size(bounds,2);
    
    out_indices = repmat((0:N_indices-1)',1,N_bounds) + repmat(bounds(1,:),N_indices,1);
    
    

end