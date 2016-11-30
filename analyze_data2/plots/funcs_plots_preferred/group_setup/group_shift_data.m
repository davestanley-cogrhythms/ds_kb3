
function gout = group_shift_data(group,shift)
    % Concatenate all the entries in group array into one large group
    % This is basically a more advanced version of group_OR; handles
    % legends better. Consider replacing?
    
    for i = 1:length(group)
        group(i).data = group(i).data + shift;
        group(i).datastats = group(i).datastats + shift;
    end
    
    gout = group;
    
end
