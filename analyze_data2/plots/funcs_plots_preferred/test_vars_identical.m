

function all_identical = test_vars_identical(varargin)
    % Tests if all input variables are identical

    var_prev = varargin{1};
    all_identical = 1;
    for i = 2:length(varargin)
        var_curr = varargin{i};
        if any(var_prev(:) ~= var_curr(:)); all_identical = 0; end
    end
    

end