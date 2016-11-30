function [group] = remove_bads_cells(bad_any,group)
    
    for i = 1:length(group)
        
        % Remove bad data
        if ~isempty(group(i).data)
            bad_data = get_grouped_data(ones(size(group(i).cells,2),1),group(i).cells,bad_any); % Pull out components of bad_any corresponding to group(i).cells
            group(i).data = group(i).data(:,~bad_data);                                         % Remove from data
        end
        
        % Remove bad stats
        if ~isempty(group(i).datastats)
            bad_data = get_grouped_data(ones(size(group(i).cells,2),1),group(i).cells,bad_any); % Pull out components of bad_any corresponding to group(i).cells
            group(i).datastats = group(i).datastats(:,~bad_data);                               % Remove from datastats
        end
        
        % Remove bad cells
        group(i).cells(bad_any,:) = false;
        
        % Redo cell counts
        group(i).numcells = sum(group(i).cells);
    end    
end
