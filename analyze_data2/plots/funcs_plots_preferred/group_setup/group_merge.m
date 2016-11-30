
function gout = group_merge(varargin)
    % Concatenate all the entries in group array into one large group
    % This is basically a more advanced version of group_OR; handles
    % legends better. Consider replacing?
    
    gr_arr = horzcat(varargin{:});
    fieldn = fieldnames(gr_arr);
    
    gout = gr_arr(1);
    for j = 1:length(fieldn)
        
%         if ~strcmp(fieldn{j},'legend') && ~strcmp(fieldn{j},'interpreter') && ~strcmp(fieldn{j},'data_name');         % Don't concatenate legend entries!
%             gout.(fieldn{j}) = cat(1,gr_arr(:).(fieldn{j}));
%             %gout.(fieldn{j}) = cat(1,gr_arr(1).(fieldn{j}));
%         end
        
        if strcmp(fieldn{j},'criteria') || strcmp(fieldn{j},'criteria_alt') || strcmp(fieldn{j},'criteria_sfc') || strcmp(fieldn{j},'ctgs') || strcmp(fieldn{j},'numcells')
            gout.(fieldn{j}) = cat(1,gr_arr(:).(fieldn{j}));
            %gout.(fieldn{j}) = cat(1,gr_arr(1).(fieldn{j}));
        end
        
    end
    
    legends_in = {gr_arr.legend};
    temp = repmat({' | '},1,length(legends_in));
    legends_in = { legends_in{:};  temp{:}};
    legends_in = legends_in(1:end-1);
    gout.legend = strcat(legends_in{:});
    
    % Shorten legend if too long
    max_legend_length = 12;
    if length(gout.legend) > max_legend_length
        gout.legend = [gout.legend(1:max_legend_length-3) '...'];
    end
    
end
