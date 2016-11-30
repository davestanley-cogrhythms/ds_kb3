

function out = is_subset(A,B)
% REturns 1 if B is a subset of A. 
% A, B are arrays of logicals (0 or 1). 
% out is a logical (0 or 1)

    if sum(A == (A | B)) == length(A)
        out = logical(1);
    else
        out = logical(0);
    end

end