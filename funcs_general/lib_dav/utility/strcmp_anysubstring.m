function out = strcmp_anysubstring(str_cells,search_string)
    % Searches through all strings str_cells; returns 1 for all cells
    % containing search_string as a substring.

    out = ~cellfun(@isempty,strfind(str_cells,search_string));

end