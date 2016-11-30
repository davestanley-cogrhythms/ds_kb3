

function group = query_legend(group,group0)
    % Send in gorup0 (usually produced by get_group_overrides2), containing
    % our desired legend entries. This function will match the ctgs in
    % group0 to the ctgs in group, and update legend entries accordignly.
    % Probably this is not right!! Should have 2 separate legends - one for
    % items grouped by criteria and one for that which is being plotted
    % (determined by ctgs).
    
    allctgs0 = [group0.ctgs];
    
    for i = 1:length(group)
        ind = allctgs0 == group(i).ctgs;
        
        if any(ind)
            group(i).legend = group0(ind).legend;
        end    
    end
    
    legends_out = [group.legend];
    
end
