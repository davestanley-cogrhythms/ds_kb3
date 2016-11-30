
function [Cave, phi,S12ave,S1ave,S2ave,F, T,confC,phistd] = coherence_spect_wrapper(varargin)

%     [P, f] = spect_wrapper(X,fs,'mode',modenum,'params',chronuxparams)
%     Calculates spectrogram using one of several
%     different methods, thus serving as a wrapper. 
%     FORMS
%         [P, f] = spect_wrapper(X)
%         [P, f] = spect_wrapper(X,options)
%     INPUTS
%         X - data (vector/matrix)
%         options - specifies options in the form of name and value pairs
%             (i.e. 'fs',1024)
%     OPTIONS
%         fs - sampling rate (Default=1)
%         mode - 1   - Default Matlab spectrogram
%                2   - Chronux (Default)
%                2.1 - Chronux Mikio (spike rate adjusted)
%         params - paramater structure to be fed into Chronux commands
%                       - note - in all cases, fs will override params.Fs.
%                       - only used if mode=2
%         Nwind - Number of datapoints to use in spectrogram window
%         fract_overlap - fraction of overlapt to use in spectrogram window
%         params - an optional structure that is passed as params to 
%                 the chronux mtspecgramc
%     OUTPUTS
%         [P, f] - power and frequency
% 
%     CONTACT: David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
% 

    % Pull inputs to local workspace
    p = inputParser;
    addRequired(p,'X',@isnumeric);
    addRequired(p,'Y',@isnumeric);
    addOptional(p,'fs',[],@isnumeric);
    addOptional(p,'mode',[],@isnumeric);     % Default mode is chronux
    addOptional(p,'Nwind',[],@isnumeric);               %  Lenght of window
    addOptional(p,'fract_overlap',[],@isnumeric);     % Fraction of overlapping amongst sliding windows
    addOptional(p,'params',[]);
    addOptional(p,'fhandle_coherence',@coherencycpb);
    parse(p,varargin{:});
    vars_pull(p.Results);
    mymode=p.Results.mode;  % Don't use mode variable name b/c overlaps mode function
    

    % Setup X
    if isvector(X); X=X(:); end
    if isvector(Y); Y=Y(:); end
    Nlfp = size(X,1);
    
    
    % Set defaults
    if isempty(fs); fs = 1; end
    if isempty(mymode); mymode = 2; end
    if isempty(Nwind); Nwind = min(round(Nlfp/10),1000); end
    if isempty(fract_overlap); fract_overlap = 0.95; end
    

    % Default params
    if isempty(params); params.tapers = [3 5];end
    
    % Parameters
    plot_on = 0;
    params.Fs = fs;
    params.err = [1 0.5];
    params.trialave = 1;
    
    
    % Check if Chronux is loaded
    if floor(mymode) == 2 && ~exist('mtspectrumc','file');
        fprintf([' Error. Chronux does not appear to be loaded. Try switching to mode = 1 to use \n' ...
            ' default Matlab commands. Alternatively, make sure that Chronux \n' ...
            ' is installed and is in your Matlab path. Exiting. \n']); return; end
    
    [~, mode_subgroups] = decimate_mode(mymode);
    mode_group = floor(mymode);
    [do_mikio] = decode_subgroups(mode_subgroups);
    

    switch mode_group;
        case 1
            Other coherency code

        case 2
            %params.tapers = get_tapers;
            %movingwind = [1 0.2];   
            Nshift = round(Nwind*(1-fract_overlap));
            
            Xstarts = 1:Nshift:(Nlfp-Nwind+1);
            Xstops = Xstarts + Nwind - 1;
            
            clear Cave phi S12ave S1ave S2ave f
            for i = 1:length(Xstarts)
                if ~do_mikio
                    [Cave{i},phi{i},S12ave{i},S1ave{i},S2ave{i},F,zerosp,confC{i},phistd{i}]=coherencycpb(X(Xstarts(i):Xstops(i),:),Y(Xstarts(i):Xstops(i),:),params);
                else
                    Ytemp = Y(Xstarts(i):Xstops(i),:);
                    StarRate = get_rateStar;
                    rate = mean(mean(Ytemp(:)))*fs;
                    [Cave{i},phi{i},S12ave{i},S1ave{i},S2ave{i},F,zerosp,confC{i},phistd{i}]=Adjcoherencycpb(X(Xstarts(i):Xstops(i),:),Ytemp,params,[],[],StarRate,rate);
                end
            end
            
            Cave = cell2mat(Cave);
            phi = cell2mat(phi);
            S12ave = cell2mat(S12ave);
            S1ave = cell2mat(S1ave);
            S2ave = cell2mat(S2ave);
            confC = cell2mat(confC);
            phistd = cell2mat(phistd);
            
            t = 1:Nlfp; t = t/fs;
            T = t(round(mean([Xstarts;Xstops])));
            
    end
    
    
    if plot_on
        %figure; imagesc([min(T),max(T)],[min(F),max(F)],log(S));set(gca,'YDir','normal')
        figure; imagesc([min(T),max(T)],[min(F),max(F)],(Cave));set(gca,'YDir','normal')
    end
end





function [mode_group, mode_subgroup] = decimate_mode(sfc_mode)  % Returns pre-decimal and post-decimal values in sfc_mode (digits as array elements)
    myprecision = 12; % If go over size of double, we'll get random numbers
    y = num2str(sfc_mode,myprecision);
    z = arrayfunu(@(x) str2num(x),y);
    dec_ind = strfind(y,'.');
    mode_group = cell2mat(z(1:dec_ind-1));
    mode_subgroup_temp = cell2mat(z(dec_ind:end));
    
    mode_subgroup = zeros(1,myprecision);
    mode_subgroup(1,1:length(mode_subgroup_temp)) = mode_subgroup_temp;
    
end



function [do_mikio] = decode_subgroups(mode_subgroups)

    if mode_subgroups(1) == 1; do_mikio = 1;
    else do_mikio = 0; end

end


