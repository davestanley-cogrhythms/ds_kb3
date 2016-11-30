
function harr = plot_scattergroupdata (x, group,sensname)

    Ngroups = length(group);
    sf = get_scaling_factor;
    
    for i = 1:length(group)
        currcell_mat = group(i).cells;
        group_cells = [];
        
        for j = 1:size(currcell_mat,2)
            group_cells = [group_cells(:); find(currcell_mat(:,j))];    % Need to make sure that the group_cells indices are ordered correctly
        end
        
        
        if sum(group_cells) > 0
            if i<=Ngroups    % Use a special symbol for all groups but the last one
                harr(i) = plot(x(group_cells),group(i).datastats,[get_clist(i) 'x'],'LineWidth',3*sf,'MarkerSize',15*sf); hold on;
                %harr(i) = plot(x(group_cells),group(i).datastats,[get_clist(i) 'o'],'LineWidth',2,'MarkerSize',17-1*i); hold on;
            else            % Use the basic symbol for the last group. 
                harr(i) = plot(x(group_cells),group(i).datastats,'k.'); hold on;
            end
        else
            harr(i) = plot(nan,[get_clist(i) 'o'],'MarkerSize',20); hold on;
        end
        
        
    end
    xlabel([sensname]); ylabel('SFC')
    
    legendarr = get_legendarr(group);
    
    legend(harr, legendarr);
    
end

