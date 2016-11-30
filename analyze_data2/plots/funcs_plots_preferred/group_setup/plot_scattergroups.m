
function harr = plot_scattergroups (x, y,group, xy_subgroup)

    if ~exist('xy_subgroup','var'); xy_subgroup = []; end
    sf = get_scaling_factor;
    
    %cells_any_group = false(size(group(1).cells,1),1);
    for i = 1:length(group)
        
        group_cells = group(i).cells;
        if isvector(group_cells); group_cells=group_cells(:);end
        group_cells = any(group_cells,2);
        
        if size(x,1) < length(group_cells)      % Hack to make code work when bad cells already removed from x and y
            if isempty(xy_subgroup); fprintf('Warning: Size of input data x/y is less than size of group.cells, but xy_subgroups is empty. \n'); end
            group_cells = group_cells(xy_subgroup);
        end
        
        if sum(group_cells) > 0
            %harr(i) = plot(x(group_cells),y(group_cells),[get_clist(i) 'o'],'LineWidth',2,'MarkerSize',16-2*i); hold on;
            %harr(i) = plot(x(group_cells),y(group_cells),[get_clist(i) 'x'],'LineWidth',14*sf,'MarkerSize',4*sf); hold on;
            %harr(i) = plot(x(group_cells),y(group_cells),[get_clist(i) 'x'],'LineWidth',1.5*sf,'MarkerSize',3*sf); hold on;
            harr(i) = plot(x(group_cells),y(group_cells),[get_clist(i) '.'],'LineWidth',14*sf,'MarkerSize',14*sf); hold on;
        else
            harr(i) = plot(0,nan,[get_clist(i) 'o'],'MarkerSize',20*sf); hold on;
        end
        
        %cells_any_group = cells_any_group | group_cells;
        
    end
    xlabel('Sens stat 1'); ylabel('Sens Stat 2')
    
    legendarr = get_legendarr(group);
    
    legend(harr, legendarr);
    
    % mmax = max([x(cells_any_group);y(cells_any_group)]);
%     plot([0 mmax],[0 mmax],'k');
    
end

