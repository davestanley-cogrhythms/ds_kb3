

function out = zscale(varargin)
    % Function that performs the compuound operation
    % z = zscale(x,args,scale_val);
    % y=zscore(x,args); z = y*scale_val;
    
    y = zscore(varargin{1:end-1});
    out = y.*varargin{end};

end