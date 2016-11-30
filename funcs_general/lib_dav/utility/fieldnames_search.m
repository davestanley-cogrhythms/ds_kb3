function out = fieldnames_search(s,search_string)
    % Searches all the field names in structure s to see if they contain
    % the text string. If they do, return only these fieldnames.

    temp = fieldnames(s);
    out = temp((~cellfun(@isempty,strfind(temp,search_string))));

end