function ismonkeyL = get_ismonkeyL(funames1D)

    % Cell dates
    if ~get_iscromer
        ismonkeyL=~cellfun(@isempty,(strfind(funames1D,'L')));
    else
        ismonkeyL=~cellfun(@isempty,(strfind(funames1D,'lu'))); % Or something like that...
    end
end