

function [S, F, T] = spect_wrapper(varargin)
%     [P, F, T] = spect_wrapper(X,fs,'mode',modenum,'params',chronuxparams)
%     Calculates spectrogram using one of several
%     different methods, thus serving as a wrapper. 
%     FORMS
%         [P, F, T] = spect_wrapper(X)
%         [P, F, T] = spect_wrapper(X,options)
%     INPUTS
%         X - data (vector/matrix)
%         options - specifies options in the form of name and value pairs
%             (i.e. 'fs',1024)
%     OPTIONS
%         fs - sampling rate (Default=1)
%         mode - 1-Default Matlab spectrogram
%                    2-Chronux (Default)
%         params - paramater structure to be fed into Chronux commands
%                       - note - in all cases, fs will override params.Fs.
%                       - only used if mode=2
%         Nwind - Number of datapoints to use in spectrogram window
%         fract_overlap - fraction of overlapt to use in spectrogram window
%         trialave - X is MxN, takes the spectrum of N columns and
%                     averages when trialave is set to 1 (defualt). Otherwise
%                     returns individual spectra for each column.
%         params - an optional structure that is passed as params to 
%                 the chronux mtspecgramc
%     OUTPUTS
%         [P, F, T] - power and frequency and time
% 
%     CONTACT: David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
% 

    % Pull inputs to local workspace
    p = inputParser;
    addRequired(p,'X',@isnumeric);
    addOptional(p,'fs',[],@isnumeric);
    addOptional(p,'mode',[],@isnumeric);     % Default mode is chronux
    addOptional(p,'trialave',[],@isnumeric);
    addOptional(p,'Nwind',[],@isnumeric);               %  Lenght of window
    addOptional(p,'fract_overlap',[],@isnumeric);     % Fraction of overlapping amongst sliding windows
    addOptional(p,'params',[]);    
    parse(p,varargin{:});
    vars_pull(p.Results);
    mymode=p.Results.mode;  % Don't use mode variable name b/c overlaps mode function
    

    % Setup X
    if isvector(X); X=X(:); end
    Nlfp = size(X,1);
    
    
    % Set defaults
    if isempty(fs); fs = 1; end
    if isempty(mymode); mymode = 2; end
    if isempty(trialave); trialave = 1; end
    if isempty(Nwind); Nwind = min(round(Nlfp/10),1000); end
    if isempty(fract_overlap); fract_overlap = 0.95; end
    

    
    % Parameters
    plot_on = 0;
    params.Fs = fs;
    params.err = [1 0.5];
    params.trialave = trialave;
    
    
    % Check if Chronux is loaded
    if mymode == 2 && ~exist('mtspectrumc','file');
        fprintf([' Error. Chronux does not appear to be loaded. Try switching to mode = 1 to use \n' ...
            ' default Matlab commands. Alternatively, make sure that Chronux \n' ...
            ' is installed and is in your Matlab path. Exiting. \n']); return; end
    
    

    switch mymode;
        case 1
            Noverlap = round(Nwind*fract_overlap);
            Ncols = size(X,2);
            S = [];
            for j = 1:Ncols
                [Stemp, F, T] = spectrogram(X(:,j),Nwind,Noverlap,Nwind,fs);
                S = cat(3,S,Stemp);
            end
                    
            if trialave; S = mean(S,3); end
            S = abs(S).^2;
            clear Stemp;

        case 2
            %params.tapers = get_tapers;
            %movingwind = [1 0.2];   
            Nshift = round(Nwind*(1-fract_overlap));
            movingwind = ([Nwind Nshift]/fs); % In units of actual time! (hopefully seconds)
            [S,T,F,Serr]=mtspecgramc(X,movingwind,params);
            S = permute(S,[2 1 3]);     % Take transpose of 1st two dimensions
    end
    
    
    if plot_on
        %figure; imagesc([min(T),max(T)],[min(F),max(F)],log(S));set(gca,'YDir','normal')
        figure; imagesc([min(T),max(T)],[min(F),max(F)],(S));set(gca,'YDir','normal')
    end
end
