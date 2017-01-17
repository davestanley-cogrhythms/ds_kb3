function out = strcmp_substr(mystr,substr)
    % Searches through all strings str_cells; returns 1 for all cells
    % containing search_string as a substring.

    if iscell(mystr) && ~iscell(substr)
        out = ~cellfun(@isempty,strfind(mystr,substr));
    elseif ~iscell(mystr) && iscell(substr)
        y = @(x) strfind(mystr,x);
        z = @(x) ~isempty(y(x));
        out = cellfun(z,substr);
    else
        out = ~isempty(strfind(mystr,search_str));
    end

end