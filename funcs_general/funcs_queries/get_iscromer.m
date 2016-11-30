


function [iscromer] = get_iscromer
    % Estimates whether we're looking at cromer data
    % based on (1) path and (2) size of sfc
    
    
    currpath = getpath('path_roy');
    iscromer = ~isempty(strfind(currpath,'cromer'));
    
    
    
end