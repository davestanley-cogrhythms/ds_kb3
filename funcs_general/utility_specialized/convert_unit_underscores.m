

function out = convert_unit_underscores (in)

    if isa(in,'str');
        out = replace_underscores(in);  % Input is a string:  in = string
    elseif iscellstr(in)
        out = cellfun(@replace_underscores,in,'UniformOutput',0);   % Input is a cell array of strings: in{1} = string
    elseif iscell(in) && iscellstr(in{1})
        out = cellfun(@(s) cellfun(@replace_underscores,s,'UniformOutput',0),in,'UniformOutput',0); %Input is of form: in{1}{1} = string
    else
        fprintf('Unknown input format. Should be either string, cell array of strings, or cell array of cell arrays of strings. \n');
    end
    
end


function out = replace_underscores(in)
    out = strrep(in,'_',' ');
end

