

function gr_out = group_OR (varargin)
    % Performs "or" operation on group criteria
    % Possibly depreciated. See group_merge

    gr_arr = horzcat(varargin{:});
    fieldn = fieldnames(gr_arr);
    
    gr_out = gr_arr(1);
    for j = 1:length(fieldn)
        gr_out.(fieldn{j}) = cat(1,gr_arr(:).(fieldn{j}));
    end
end