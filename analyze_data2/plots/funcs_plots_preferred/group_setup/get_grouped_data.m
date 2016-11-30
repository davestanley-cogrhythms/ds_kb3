
function [group_data, group_funames1D] = get_grouped_data(ctgs,group_cells,data,funames)
    
    statfun = get_statfunc;
    
    plot_on = 0;


    Nsubgroups = size(ctgs,1);
    Nsubgroups2 = size(group_cells,2);
    
    
    if Nsubgroups ~= Nsubgroups2
        fprintf('Error - rows in ctgs should match number of columns in group_cells. Ending \n');
        group_data = nan;
        return;
    end
    
    Ndata = size(data,1);
    cells_per_subgroup = sum(group_cells);
    cells_cum = cumsum(cells_per_subgroup);
    cells_cum = [0 cells_cum];
    
    tot_cells = sum(cells_per_subgroup);
    N4D = size(data,4);
    group_data = zeros(Ndata,tot_cells,N4D)*nan;
    
    
    for i = 1:Nsubgroups
        if any(ctgs(i,:) > size(data,3)); warning('Current ctgs out of data index. Probably means that group.ctgs is too large.'); end
        curr_cells = group_cells(:,i);
        group_data(:,cells_cum(i)+1:cells_cum(i+1),:) = squeeze(statfun(data(:,curr_cells,ctgs(i,:),:),3));
        if nargout > 1; group_funames{i} = funames(curr_cells); end
    end
    
    if nargout > 1; group_funames1D = cat(2,group_funames{:}); end
    
end