

function [h,hl] = plott_all(varargin)
%     [h, hl] = plott_all(varargin)
%     Plots the time series, psd, and spectrogram in a series of subplots.
%     FORMS
%         [h, hl] = plott_all(X)
%         [h, hl] = plott_all(t,X)
%         [h, hl] = plott_all(X,options)
%         [h, hl] = plott_all(t,X,options)
%     INTPUTS
%         t - 1D vector of times.
%         X - data (vector/matrix)
%     OPTIONS
%         fs - sampling rate (Default=1; when t is supplied fs is overridden by spacing of t )
%         mode - 1-Default Matlab spectra commands (Pwelch, etc)
%                2-Chronux (Default)
%         params - paramater structure to be fed into Chronux commands
%                       - note - in all cases, fs will override params.Fs.
%                       - only used if mode=2
%         logplot - log scale plotting (boolean, Default=0)
%         psd_on - Set to 1/0 to display/hide psd. This can also be
%               a cell array of name value option pairs to pass to plott_psd
%         spect_on - Set to 1/0 to display/hide spectrogram. This can also be
%               a cell array of name value option pairs to pass to plot_spec 
%         Nwind - Number of datapoints to use in spectrogram window
%         axis_lims - optional 2 element matrix that specifies range of frequency axis
%                 (i.e. [0 200] = 0 to 200 Hz). Default is [].
%         fsubplot - subplot function to use. Default is @subplotsq
%     OUTPUTS
%         h = collection of subplot object handles
%         hl = line handles from each subplot
% 
%     Example 1 - Basic plot
%     load pinknoise
%     X = x(1:1000); t = 1:size(X,1);
%     figure; out2 = plott_all(X,'fs',1e3,'mode',1,'logplot',1);
%     
%     Example 2 - Plot of 2D matrix
%     X = cumsum(randn([1000,5]));
%     t = 1:size(X,1);
%     figure; out2 = plott_all(X,'fs',1e3,'mode',1,'logplot',1);
%     
%     Example 3 - Passing options to spect on
%     load chirp
%     Fs = 500;
%     t = 1:length(y); t=t/Fs;
%     figure; plott_all(y,'fs',Fs,'mode',1,'spect_on',{'Nwind',300},'logplot',1);
% 
%     CONTACT: David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
% 



    % Pull inputs to local workspace
    [reg, args]=parseparams(varargin);  % Reg should contain t and X
    p = inputParser;
    addOptional(p,'fs',1,@isnumeric);
    addOptional(p,'mode',2,@isnumeric);     % Default mode is chronux
    addOptional(p,'logplot',0,@isnumeric);
    addOptional(p,'psd_on',1);
    addOptional(p,'spect_on',1);
    addOptional(p,'fsubplot',@subplotsq); 
    addOptional(p,'params',[]); 
    addOptional(p,'axis_lims',[],@isnumeric);
    addOptional(p,'Nwind',[],@isnumeric);
    parse(p,args{:});
    vars_pull(p.Results)
    psd_mode = p.Results.mode;
    
    % Set defaults
    if isempty(fs); fs = 1; end
    
    % Load X
    [t,X,fs] = sort_out_timeseries(reg,fs); % Assigns values to t, X, and fs. If available, the fs from t always overrides fs.
    
    % Set arguments for passing to spect and psd functions.
    psd_args={};spect_args={};
    if iscell(psd_on); psd_args = psd_on; psd_on = 1; end
    if iscell(spect_on); spect_args = spect_on; spect_on = 1; end
    
    cp=1; % Current plot
    fname{cp}=@(varargin) plot_with_labels(varargin{:});
    
    
    if spect_on
        cp=cp+1;
        fname{cp}=@(t,x) plott_spect_cust(t,x,'fs',fs,'mode',mode,'logplot',logplot,'axis_lims',axis_lims,'Nwind',Nwind,'params',params,spect_args{:});
    end
    
    if psd_on
        cp=cp+1;
        fname{cp}=@(t,x) plott_psd_cust(t,x,'fs',fs,'mode',psd_mode,'logplot',logplot,'axis_lims',axis_lims,'params',params,psd_args{:});
    end
    
    [h,hl] = plott_handles(fname,fsubplot,t,X);
    
%     
%     
%     h(cp) = subplotsq(Nplots,cp);
%     hl{cp} = plot(t,X);
%     xlabel('Time');
%     
%     if plotpsd_on
%         cp=cp+1;
%         h(cp)=subplotsq(Nplots,cp);
%         
%         hl{cp} = plott_psd(X,1/dt,psd_mode);
% 
%     end
%     
%     if plotspect_on
%         cp=cp+1;
%         [P, f, t] = spect_wrapper(X(:,1),1/dt, psd_mode);
%         h(cp) = subplotsq(Nplots,cp);
%         
%         if ~logplot
%             hl{cp} = imagesc([min(T),max(T)],[min(F),max(F)],(S));set(gca,'YDir','normal')
%         else
%             hl{cp} = imagesc([min(T),max(T)],[min(F),max(F)],log(S));set(gca,'YDir','normal')
%         end
%         
%         xlabel('Time');
%         ylabel('Freq');
%         colorbar
%     end

end


% 
% function h=subplotsq(N,i)
% 
%     Ncols = ceil(sqrt(N));
%     Nrows = ceil(N/Ncols);
%     h=subplot(Nrows,Ncols,i);
%     
% end
% 
% 
% 
% 
% function tapers = get_tapers
%     % Returns the default value of tapers that are used throughout
%     % the code.
%     tapers = [3 5];
%     %tapers = [5 9];
% end





function h = plot_with_labels(varargin)

    t=varargin{1}; x=varargin{2};
    
%     plot(varargin{:});
        
    opt.shift = 0.0;
    h = plot_matrix2(t,zscore(x),opt);
    
    if isvector(t); t=t(:); end
    tmean = mean(t,2);
    
    xlabel('time');
    xlim([min(tmean) max(tmean)]);
    title('Data');
    xlabel('time (s)');

end


function h = plott_psd_cust(varargin)

    h = plott_psd(varargin{:});
    %xlim([0 500]);
    title('PSD');
    xlabel('Freq (Hz)');
    ylabel('PSD');

end



function h = plott_spect_cust(varargin)

    h = plott_spect(varargin{:});
    %ylim([0 500]);
    title('Spect');
    ylabel('Freq (Hz)');
    xlabel('Time(s)');

end