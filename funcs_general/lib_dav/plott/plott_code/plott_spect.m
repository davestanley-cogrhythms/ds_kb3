

function hl = plott_spect(varargin)
%     [hl] = plott_spect(X,fs,'mode',modenum,'params',chronuxparams,'logplot',logplot)
%     Takes input time series and plots the spectrogram using one of several methods
%     (calls spect_wrapper.m)
%     v1.0 David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
%     FORMS
%         [hl] = plott_spect(X)
%         [hl] = plott_spect(t,X)
%         [hl] = plott_spect(X,options)
%         [hl] = plott_spect(t,X,options)
%     INPUTS
%         t - times (vector)
%         X - data (vector/matrix)
%         options - specifies options in the form of name and value pairs
%             (i.e. 'fs',1024)
%     OPTIONS
%         fs - sampling rate (Default=1; when t is supplied fs is overridden by spacing of t )
%         mode - 1-Matlab spectrogram function
%                2-Chronux (default, must have Chronux in path)
%                3-Chronux (old implementation, using spect_wrapper function)
%         params - paramater structure to be fed into Chronux commands
%                       - note - in all cases, fs will override params.Fs.
%                       - only used if mode=2
%         Nwind - Number of datapoints to use in spectrogram window
%         fract_overlap - fraction of overlapt to use in spectrogram window
%         logplot - log scale plotting (boolean, Default=0)
%         axis_lims - optional 2 element matrix that specifies range of frequency axis
%                 (i.e. [0 200] = 0 to 200 Hz). Default is [].
%         clims - optional 2 element matrix that specifies range of colors
%                 (i.e. [0 200] = 0 to 200 Hz). Default is [].
%         zscoreplot - takes zscore along rows of matrix, effectively
%                     normalizing the power in all frequency bands (boolean, Default=0)
%     OUTPUTS
%         [HL] - plot handle
% 
%     
%     Example 1 - Basic plot
%     load pinknoise
%     X = x;
%     t = 1:size(X,1);
%     subplotrows(2,1); out1 = plot(t,X);
%     subplotrows(2,2); out2 = plott_spect(X,'fs',1e3,'mode',1,'logplot',1);
% 
%     Example 2 - Plot 2D Matrix
%     X = cumsum(randn([1000,5]));
%     t = 1:size(X,1);
%     subplotrows(2,1); out1 = plot(t,X);
%     subplotrows(2,2); out2 = plott_spect(X,'fs',1e3,'mode',1,'logplot',1);
% 
%     Example 3
%     fs = 1e3;
%     t=(0:1/fs:1)'; X = [sin(2*20*pi*t.*t)]+0;
%     subplotrows(2,1); out1 = plot(t,X);
%     subplotrows(2,2); out2 = plott_spect(X,'fs',fs,'mode',1,'Nwind',300); ylim([0 100])
% 
%     Example 4
%     load chirp
%     t = 1:length(y); t=t/Fs;
%     subplotrows(2,1); out1 = plot(t,y);
%     subplotrows(2,2); out2 = plott_spect(y,'fs',Fs,'mode',2,'Nwind',300,'logplot',1);
% 
%     CONTACT: David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
% 

    % Pull inputs to local workspace
    [reg, args]=parseparams(varargin);  % Regs should contain t (optional) and X
    p = inputParser;
    addOptional(p,'fs',1,@isnumeric);
    addOptional(p,'mode',[],@isnumeric);     % Default mode is chronux
    addOptional(p,'logplot',[],@isnumeric);
    addOptional(p,'zscoreplot',[],@isnumeric);
    addOptional(p,'Nwind',[],@isnumeric);               %  Lenght of window
    addOptional(p,'fract_overlap',[],@isnumeric);       % Fraction of overlapping amongst sliding windows
    addOptional(p,'params',[]); 
    addOptional(p,'axis_lims',[],@isnumeric);
    addOptional(p,'clims',[],@isnumeric);
    parse(p,args{:});
    vars_pull(p.Results);
    psd_mode = p.Results.mode;
    
    %Set defaults
    if isempty(fs); fs = 1; end
    if isempty(logplot); logplot = 0; end
    if isempty(logplot); zscoreplot = 0; end
    if isempty(psd_mode); psd_mode = 2; end
    
    %Setup X
    [t,X,fs] = sort_out_timeseries(reg,fs); % Assigns values to t, X, and fs. If available, the fs from t always overrides fs.
    t = mean(t,2); % Convert t back to vector - we don't need it to be a matrix in this case
    
    if psd_mode == 1
        [S, F, T] = spect_wrapper(X,fs,'mode',1,'params',params,'trialave',1,'Nwind',Nwind,'fract_overlap',fract_overlap);
    elseif psd_mode == 2
        if isempty(params);
            params.tapers = [3,5];
            params.Fs = fs;
        end
        fname = @(X) mtspectrumc(X,params);
        [S, F, T] = func2spectrogram(X,'fname',fname,'fs',fs,'trialave',1,'Nwind',Nwind,'fract_overlap',fract_overlap);     % IUse the more general function, based on function handles
    elseif psd_mode == 3
        [S, F, T] = spect_wrapper(X,fs,'mode',2,'params',params,'trialave',1,'Nwind',Nwind,'fract_overlap',fract_overlap);   % Old way of doing things.
        
    end
        
    
    
    
    
    
    T = T - mean(T) + mean(t) + 1/fs/2; % Shifts the spectrogram to be centered on wherever the original input t was centred.
    
    if zscoreplot
        S = zscore(S,[],2);
    end

    if ~logplot
        hl = imagesc([min(T),max(T)],[min(F),max(F)],(S));set(gca,'YDir','normal');
    else
        hl = imagesc([min(T),max(T)],[min(F),max(F)],10*log10(S));set(gca,'YDir','normal');
    end

    if ~isempty(clims); caxis(clims);end
    xlabel('Time');
    ylabel('Freq');
    xlim([min(t) max(t)]);
    if ~isempty(axis_lims)
        ylim(axis_lims);
    end

end





