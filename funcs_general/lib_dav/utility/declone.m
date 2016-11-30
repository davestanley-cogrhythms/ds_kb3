

function out = declone(x1,x2)
    % Ensures that x1 equals x2 and returns x1
    % If x1 ~= x2, returns error
    % If either x1 or x2 are empty, ignores them
    
    
    if ~isempty(x1) && ~isempty(x2)
        if any(x1(:) ~= x2(:)); error('Mismatch between x1 and x2'); end
        out = x1;
    elseif ~isempty(x1)
        out = x1;
    elseif ~isempty(x2)
        out = x2;
    else
        out = [];
    end
    

end