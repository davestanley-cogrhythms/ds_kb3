

function legend_arr = get_legendarr(group,detail_level,show_ncells,force_plaintext)
% legend_arr = get_legendarr(group)

    if nargin < 4
        force_plaintext = 0;
    end
    if nargin < 3
        show_ncells = 1;
    end
    if nargin < 2
        detail_level = Inf;
    end
    

    for i = 1:length(group)
        
        legend_arr{i} = ''; % Empty
        
        if detail_level >= 1
            legend_arr{i} = [num2str(i)];   % Just the number
        end
        
        if detail_level >= 2    % "Gr1"
            legend_arr{i} = ['Gr' num2str(i)];
        end
        
        if detail_level >= 3    % Inherit name from group
            
            % First import plaintext legend if it exists
            if ~isempty(group(i).legend_plaintext)
                legend_arr{i} = group(i).legend_plaintext;
            end
            
            % Override if we're not forcing plaintext
            if ~isempty(group(i).legend) && ~force_plaintext
                legend_arr{i} = group(i).legend;
            end
            
        end
        
        if show_ncells    % All of the above, plus
            %cellnums_string = num2str(group(i).numcells);
            cellnums_string = num2str(sum(group(i).numcells));
            legend_arr{i} = [legend_arr{i} ' N=' cellnums_string];   % Add on cell numbers!
        end
        
    end
    
    

end