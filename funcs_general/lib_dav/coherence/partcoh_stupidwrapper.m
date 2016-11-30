

function [varargout] = partcoh_stupidwrapper(varargin)

    check_data_identical = 0; % Some error checking, but can slow things down.

    if nargin < 3
        warning('Must have a params argument to contain fname and chosen_dims metadata.');
    end
    data1_curr = varargin{1};
    data2_curr = varargin{2};
    params = varargin{3};
    
    % Any extra argumetns we might have
    if nargin > 3
        args = varargin(4:end);
    else args = {};
    end
    
    fname = params.fname;
    chosen_dims = params.chosen_dims;
    
    if check_data_identical
        if any(data1_curr(:) ~= data2_curr(:))
            warning('For this stupid wrapper, I was expecting data1_curr == data2_curr. This is not a big deal, but should investigate.');
        end
    end
    
    [varargout{1:nargout}] = partcoh(fname,data1_curr,chosen_dims,params,args{:});
    %[varargout{1:nargout}] = partcoh_full(fname,data1_curr,params,args{:});


end

