

function hl = plott_psd(varargin)
%     [hl] = plott_psd(X,fs,'mode',mode,'params',params,'logplot',logplot)
%     Takes input time series and plots the PSD using one of several methods
%     (calls psd_wrapper.m)
%     FORMS
%         [hl] = plott_psd(X)
%         [hl] = plott_psd(t,X)
%         [hl] = plott_psd(X,options)
%         [hl] = plott_psd(t,X,options)
%     INPUTS
%         t - times (vector)
%         X - data (vector/matrix)
%         options - specifies options in the form of name and value pairs
%             (i.e. 'fs',1024)
%     OPTIONS
%         fs - sampling rate (Default=1; when t is supplied fs is overridden by spacing of t )
%         mode - 1-Matlab pWelch
%                2-Chronux (Default, must have Chronux in your path)
%         params - paramater structure to be fed into Chronux commands
%                       - note - in all cases, fs will override params.Fs.
%                       - only used if mode=2
%         logplot - log scale plotting (boolean, Default=0)
%         axis_lims - optional 2 element matrix that specifies range of frequency axis
%                 (i.e. [0 200] = 0 to 200 Hz). Default is [].
%     OUTPUTS
%         [HL] - plot handle
% 
%     Example 1 - Baisc plot
%     load pinknoise
%     X = x;
%     t = 1:size(X,1);
%     subplotrows(2,1); out1 = plot(t,X);
%     subplotrows(2,2); out2 = plott_psd(X,'fs',1e3,'mode',1,'logplot',1);
% 
%     Example 2 - Plot 2D matrix
%     X = cumsum(randn([1000,5]));
%     out = plott_psd(X,'fs',1e3,'mode',1,'logplot',1);
%     
%     Example 3 - Plot 2D matrix in Chronux mode, with passing 
%                 params structure to Chronux
%     fs = 1e3;
%     t=(0:1/fs:1)'; X = [sin(2*20*pi*t) 1.1*sin(2*21*pi*t)]+1.1;
%     params.trialave = 1;
%     subplotrows(2,1); out1 = plot(t,X);
%     subplotrows(2,2); out = plott_psd(X,'fs',fs,'mode',2,'params',params);
% 
%     CONTACT: David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
% 

    
    % Pull inputs to local workspace
    [reg, args]=parseparams(varargin);  % Regs should contain t and X
    p = inputParser;
    addOptional(p,'fs',[],@isnumeric);
    addOptional(p,'mode',[],@isnumeric);     % Default mode is chronux
    addOptional(p,'logplot',[],@isnumeric);
    addOptional(p,'params',[]);
    addOptional(p,'plotargs',{});
    addOptional(p,'axis_lims',[]);
    parse(p,args{:});
    vars_pull(p.Results);
    psd_mode = p.Results.mode;
    
    % Set Defaults
    if isempty(fs); fs = 1; end
    if isempty(logplot); logplot = 0; end
    if ~isfield(params,'trialave'); params.trialave = 1; end
    
    
    % Setup X
    [t,X,fs] = sort_out_timeseries(reg,fs); % Assigns values to t, X, and fs. If available, the fs from t always overrides fs.
    

    [P, f] = psd_wrapper(X,fs,'mode',psd_mode,'params',params);
        
    if ~logplot
        hl = plot(f,P,plotargs{:});
    else
        hl = loglog(f,P,plotargs{:});
        set(gca,'XScale','log')
        set(gca,'YScale','log')
    end
    xlabel('Freq');
    if ~isempty(axis_lims)
        xlim(axis_lims);
    end
    
end

