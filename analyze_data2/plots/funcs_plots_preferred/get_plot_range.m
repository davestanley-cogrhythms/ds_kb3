




function [gmin, gmax] = get_plot_range(datacell)
%     Searches through cells in datacell and
%     returns the max and min values across all cells

    maxes = cellfun(@(x) max(x(:)),datacell);
    mins = cellfun(@(x) min(x(:)),datacell);
    
    gmin = min(mins);
    gmax = max(maxes);
    
end


% function [gmin, gmax] = get_plot_range(group)
%     % As above, but works instead on group
% 
%     
%     all_data = {group.data};
%     ad = cellfunu(@(x) mean(x,2), all_data);
% 
%     adall = cat(2,ad{:});
%     gmin = min(adall(:));
%     gmax = max(adall(:));
% end