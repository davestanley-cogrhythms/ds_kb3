

function hl = plott_coherency(varargin)
%     [hl] = plott_coherency(X,Y,fs,'mode',modenum,'params',chronuxparams,'logplot',logplot)
%     Takes input time series and plots the spectrogram using one of several methods
%     (calls spect_wrapper.m)
%     v1.0 David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
%     FORMS
%         [hl] = plott_coherency(X,Y)
%         [hl] = plott_coherency(t,X,Y)
%         [hl] = plott_coherency(X,Y,options)
%         [hl] = plott_coherency(t,X,Y,options)
%     INPUTS
%         t - times (vector)
%         X - data (vector/matrix)
%         Y - data (vector/matrix)
%         options - specifies options in the form of name and value pairs
%             (i.e. 'fs',1024)
%     OPTIONS
%         fs - sampling rate (Default=1; when t is supplied fs is overridden by spacing of t )
%         mode - 1-Matlab spectrogram function (not implemented)
%                2-Chronux (default, must have Chronux in path)
%         params - paramater structure to be fed into Chronux commands
%                       - note - in all cases, fs will override params.Fs.
%                       - only used if mode=2
%         Nwind - Number of datapoints to use in spectrogram window
%         fract_overlap - fraction of overlapt to use in spectrogram window
%     OUTPUTS
%         [HL] - plot handle
% 
%     Example 1 - Decreasing high frequency coherence
%     % Signal X is pink noise
%     load pinknoise
%     X = x;
%     X = reshape(X,10000,10);
%     
%     % Signal Y adds some random white noise to X of ever increasing
%     % strength, decreasing coherence at high frequencies.
%     scale = repmat(linspace(0.005,0.02,length(X))',1,size(X,2));
%     Y = X+randn(size(X)).*scale;
%     t = 1:size(X,1);
%     subplotsq(3,1); out1 = plot(t,[mean(X,2) mean(Y,2)+.02]);
%     subplotsq(3,2); out2 = plott_coherency(X,Y,'fs',1e3); colorbar;
%     subplotsq(3,3);
%         params.trialave = 1;
%         [Cave,phi,S12ave,S1ave,S2ave,f]=coherencycpb(X,Y,params);
%         plot(f,Cave)
    
%     CONTACT: David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
% 

    % Pull inputs to local workspace
    [reg, args]=parseparams(varargin);  % Regs should contain t (optional) and X
    p = inputParser;
    addOptional(p,'fs',1,@isnumeric);
    addOptional(p,'mode',[],@isnumeric);     % Default mode is chronux
    addOptional(p,'Nwind',[],@isnumeric);               %  Lenght of window
    addOptional(p,'fract_overlap',[],@isnumeric);       % Fraction of overlapping amongst sliding windows
    addOptional(p,'params',[]); 
    parse(p,args{:});
    vars_pull(p.Results);
    psd_mode = p.Results.mode;
    
    %Set defaults
    if isempty(fs); fs = 1; end
    if isempty(params)
        params.Fs = fs;
        params.tapers = [3,5];
        params.trialave = 1;
    end
    
    %Setup X
    [t,X,Y,fs] = sort_out_timeseries2(reg,fs);
    t = mean(t,2); % Convert t back to vector - we don't need it to be a matrix in this case
    
    fname = @(X,Y) coherencycpb(X,Y,params);
    [Cave,~,~,~,~,F,T] = func2spectrogram(X,Y,'fs',fs,'Nwind',Nwind,'fract_overlap',fract_overlap,'fname',fname,'trialave',1);
%     [Cave,~,~,~,~,F,T] = coherence_spect_wrapper(X,Y,'fs',fs,'mode',psd_mode,'params',params,'Nwind',Nwind,'fract_overlap',fract_overlap);
    T = T - mean(T) + mean(t) + 1/fs/2; % Shifts the spectrogram to be centered on wherever the original input t was centred.
    
    

    hl = imagesc([min(T),max(T)],[min(F),max(F)],(Cave));set(gca,'YDir','normal');

    xlabel('Time');
    ylabel('Freq');
    xlim([min(t) max(t)]);
    %ylim([0 100]);

end





function [t,X,Y,fs] = sort_out_timeseries2(reg,fs)
    if length(reg) == 2
        X=reg{1}; if isvector(X); X=X(:); end
        Y=reg{2}; if isvector(Y); Y=Y(:); end
        t=1:size(X,1); t=t/fs;
    elseif length(reg) == 3
        t=reg{1}; 
        X=reg{2}; if isvector(X); X=X(:); end
        Y=reg{3}; if isvector(Y); Y=Y(:); end
        fs = 1./(mode(diff(t)));
        fs = fs(1);             % Take first entry if t is a matrix
    else
        fprintf('too many numerical inputs \n');
    end
    
    % Ensure that t is a column matrix
    if isvector(t);
        t=t(:); 
        
        % If X is a matrix, expand t to also be a matrix of the same dim
        if ~isvector(X);
            t = repmat(t,1,size(X,2));
        end
    end
end


