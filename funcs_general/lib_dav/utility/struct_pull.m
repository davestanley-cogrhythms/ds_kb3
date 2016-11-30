function [varargout] = struct_pull(s,fld)
    % This is depreceated - use vars_pull
    
    
    % If no flds specified, just take all
    if nargin < 2
        fld = fieldnames(s{1});
    end
    
    % If fld is not a cell array, make it one.
    if ~iscell(fld)
        fld = {fld};
    end

    
    for i = 1:length(fld)
        name = fld{i};
        if iscell(s)
            varargout{i} = cellfun(@(x) x.(name),s,'UniformOutput',0);
        else
            varargout{i} = arrayfun(@(x) x.(name),s,'UniformOutput',0);
        end
        
    end
    
    % If output is only 1 dimensional convert it from a cell to a matrix
%     if length(varargout) == 1
%         varargout = varargout{1};
%     end
    
end

