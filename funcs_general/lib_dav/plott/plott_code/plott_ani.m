
function [returns] = plott_ani(varargin)
%     [returns] = plotani(X,fs,options)
%     Takes in a 2D or 3D matrix and generates an animated plot that cycles through
%     the matrix columns when the user presses "enter." If multiple function
%     handles are specified in a cell array, then each function handle
%     will get its own subplot.
%     FORMS
%         [returns] = plott_ani(X)
%         [returns] = plott_ani(t,X)
%         [returns] = plott_ani(X,options)
%         [returns] = plott_ani(t,X,options)
%     INTPUTS
%         t - 1D vector of times. This is just used to calculate fs
%         X - vector or matrix of data, where rows are times and columns are variables
%           - if X is 3D, then the 2nd dimension is used for animation and the 3rd is plotted with hold on
%         options - specifies some plotting options in the form of name and value pairs
%     OPTIONS
%           'fs',fs - sampling rate (Default=1; when t is supplied fs is overridden by spacing of t )
%           'randcol', boolean -  randomize columns for plotting. BOOLEAN = 1 for true
%           'fname', fname - function handle or cell array of multiple function handles
%           'splitup',N - automatically converts large datafiles into smaller chuncks of size N
%                 by "chunking" along the 1st dimension. See examples. Zero for no
%                 splitting (default).
%           'fsubplot', subplotting function to use. Default @subplotsq
%           'plotargs' - cell array of option pairs to be passed to all plotting functions
%     OUTPUTS
%         returns - set of user keystrokes
% 
%     Example 1 - Plot a 3D matrix
%     X = abs(randn([100,5,2])); fname = {@plot, @loglog};
%     t = 1:100; t = t/100;
%     out = plott_ani(t,X,'fs',1e-3,'randcol',0,'fname',fname);
%     
%     Example 2 - Plot a 2D matrix that is converted to 3D using splitup
%     t=(1:950)'; X = [sin(2*pi/21*t) 1.1*sin(2*pi/21*t)]+1.1;
%     fname = {@plott_fs, @loglog, @(x) plott_psd(x,'fs',1)};
%     out = plott_ani(X,'fname',fname,'splitup',100);
% 
%     CONTACT: David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
%     

    [reg, args]=parseparams(varargin);  % Regs should contain t and X
    p = inputParser;
    addOptional(p,'fs',[],@isnumeric);
    addOptional(p,'randcol',0,@isnumeric);
    addOptional(p,'splitup',0,@isnumeric);
    addOptional(p,'meanMode',0,@isnumeric);
    addOptional(p,'fname',@plot);
    addOptional(p,'fsubplot',@subplotsq);
    addOptional(p,'plotargs',{});
    
    parse(p,args{:});
    vars_pull(p.Results)
    
    % Set Defaults
    if isempty(fs); fs = 1;
    else plotargs = {plotargs{:}, 'fs',fs};
    end
    

    [t,X,fs] = sort_out_timeseries(reg,fs); % Assigns values to t, X, and fs. If available, the fs from t always overrides fs.
    
    if isvector(X);
        X=X(:);
        if splitup <= 0 splitup = round(length(X)/10); end % Turn on splitup, since it's pointless otherwise
    end
    
    if splitup > 0
        X = splitup_data(X,splitup);     % Take the data in each row, and split it into columns, which will then be animated across
        if isnan(sum(X(:))); returns=NaN;return; end % splitup returns NaN if X is wrong dimensionality.
    end
    
    if meanMode
        X=mean(X,2);
    end

    N = size(X,1);
    Nepochs = size(X,2);
    
    if randcol cols = randperm(Nepochs);
    else cols = 1:Nepochs; end
    
    %t=1:N;t=t/fs;
    
%     for i = 3
    for i = 1:length(cols)
        clf
        if ndims(X) <= 2
            Xcurr = (X(:,cols(i)));
        elseif ndims(X) == 3
            Xcurr = squeeze(X(:,cols(i),:));
        end
        
        if length(reg) == 1
            plott_handles(fname,fsubplot,Xcurr,plotargs{:});
        elseif length(reg) == 2
            t = 1:size(Xcurr,1); t = t / fs;
            plott_handles(fname,fsubplot,t,Xcurr,plotargs{:});
        end
        
        prompt = 'Type a comment or hit enter to continue q to quit: ';
        returns{i} = input(prompt,'s');
%         returns{i}=1;
        
        if strcmp(returns{i},'q') || strcmp(returns{i},'Q'); break; end
        
    end
    
end


function vars_pull(s)
    for n = fieldnames(s)'
        name = n{1};
        value = s.(name);
        assignin('caller',name,value);
    end
end


