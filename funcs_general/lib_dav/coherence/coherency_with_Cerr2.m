

function [varargout] = coherency_with_Cerr2(varargin)
%     Usage: [varargout] = coherency_with_Cerr(fname,data1,data2,params,args)
%     Implements [varargout] = fname(data1,data2,params,args{:}). If there are too
%     many output arguments requested and jackknife is turned off, it then
%     generates fake data (zeros) to return for the jackknife. This is
%     useful when you want to have generic code that will run regardless of
%     whether we are in jackknife mode.
%     Forms:
%         [varargout] = coherency_with_Cerr(fname,data1,data2,params,args)
% 
%     The older version of this code is more generic, but unfortunately it
%     must run everything twice if it fails to return the appropriate input
%     the first time.


    fname = varargin{1};
    data1 = varargin{2};
    data2 = varargin{3};
    
    if nargin > 3
        params = varargin{4};
    else
        params = [];
    end
    
    if nargin > 4
        args = varargin(5:end);
    else
        args = {};
    end
    
    
    do_jackknife = 0;
    if isfield(params,'err');
        if params.err(1) == 2
            do_jackknife = 1;
        end
    end
    
    
    if nargout >= 9 && do_jackknife
        [varargout{1:nargout}] = fname(data1,data2,params,args{:});
    elseif nargout >= 9 && ~do_jackknife
        [varargout{1:nargout-1}] = fname(data1,data2,params,args{:});
        Cerr1_fake =  zeros(1,length(varargout{1}(:)));
        Cerr2_fake =  ones(1,length(varargout{1}(:)));          % Needs to match the dimensionality of Cave for cases where we do skip_Jackknife_for_ctg_alltrials17
        varargout{nargout} = [Cerr1_fake; Cerr2_fake];
    else
        [varargout{1:nargout}] = fname(data1,data2,params,args{:});
    end
    

end