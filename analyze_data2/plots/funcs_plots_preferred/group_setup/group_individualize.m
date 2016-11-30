


function gr_sp = group_individualize(group)
    % This code splits a single input group, containing multiple cells,
    % into a separate group for each cell.
    
    if length(group) > 1
        error('Can only take in a group array of size 1.');
    end
    
    % Extract cells info from group
    all_cells = find(group.cells);
    
    % Build template from which we will derive gr_out
    grt = group;
    grt.cells(:) = false;
    grt.data = [];
    grt.datastats = [];
    grt.funames = {};
    grt.numcells = 1;
    
    % Initialize array gr_out with values from template
    Nc = group.numcells;
    gr_sp = repmat(grt,1,Nc);       % Split version of group
    
    % Import values for data, datastats, and funames
    for i = 1:Nc
        gr_sp(i).cells(all_cells(i)) = true;
        gr_sp(i).data = group.data(:,i,:);
        gr_sp(i).datastats = group.datastats(i);
        gr_sp(i).funames = {group.funames{i}};
    end
    
end
