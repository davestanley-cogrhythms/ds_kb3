

function [varargout] = func2spectrogram(varargin)
%     [varargout] = func2spectrogram(X,..,OPTIONS)
%     Takes data and a function handle as input, and then
%     iteratively applies that function handle over segments
%     of the data (of size Nwind, with overlap fract_overlap).
%     The number of data arguments passed must match the number
%     of data arguments taken by the function handle. Likewise
%     the number of outputs returned equals the number of 
%     outputs from the function handle, plus 1 for the vector T
%     indicating time bins of the spectrogram.
%     FORMS
%         [fout, T] = func2spectrogram(X)
%         [fout, T] = func2spectrogram(X,Y)
%         [fout, T] = func2spectrogram(X,Y,Z,...)
%         [fout, T] = func2spectrogram(X,...,OPTIONS)
%     INPUTS
%         X       - data (vector/matrix)
%         OPTIONS - specifies options in the form of name and value pairs
%             (i.e. 'fs',1024)
%     OPTIONS (name value pairs)
%         fname - Function handle. Needs to be of the form:
%                     [fout] = fname(X)
%                     [fout] = fname(X,Y)
%                     [fout] = fname(X,Y,Z,...)
%                 where fname is some function that calcualtes spectrogram
%                 where fout is the cell array of fname outputs
%                 and ndims(fout{:}) should be <= 2.
%         fs - sampling rate (Default=1)
%         Nwind - Number of datapoints to use in spectrogram window
%         fract_overlap - fraction of overlapt to use in spectrogram window
%         trialave - X is MxN, takes the spectrum of N columns and
%                     averages when trialave is set to 1 (defualt). Otherwise
%                     returns individual spectra for each column.
%     OUTPUTS
%         fout      - outputs of fname function (i.e. spectral power), frequency, and time
%         T         - Time points of spectrogram
% 
%     CONTACT: David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
% 

    % Parameters
    plot_on = 0;
    
    % Pull inputs to local workspace
    [X, args]=parseparams(varargin);  % Regs should contain t (optional) and X
    p = inputParser;
    addOptional(p,'fname',[]);
    addOptional(p,'fs',[],@isnumeric);
    addOptional(p,'trialave',[],@isnumeric);
    addOptional(p,'Nwind',[],@isnumeric);               %  Lenght of window
    addOptional(p,'fract_overlap',[],@isnumeric);     % Fraction of overlapping amongst sliding windows
    parse(p,args{:});
    vars_pull(p.Results);
    

    % Setup X
    for i = 1:length(X);
        if isvector(X{i}); X{i} = X{i}(:); end
    end
    Nlfp = size(X{1},1);
    
    % Set defaults
    nout = nargout;
    indexF = 2;
    if isempty(fname);
        params.Fs = fs;
        params.tapers = [3,5];
        params.trialave = 0;
        if length(X) == 2
            fname = @(X,Y) coherencycpb(X,Y,params); % Autoguesses spike-field coherence if 2 inputs given
            indexF = 6;                              % Index of output argument holding F. This is only used for plot_on
            nout = max(indexF+1,nargout);
        elseif length(X) == 1
            fname = @(X) mtspectrumc(X,params); % Autoguesses spectrogram if 1 input given
            indexF = 2;
            nout = max(indexF+1,nargout);
        end
    end
    
    if isempty(fs); fs = 1; end
    if isempty(Nwind); Nwind = min(round(Nlfp/10),1000); end
    if isempty(fract_overlap); fract_overlap = 0.9; end
    if isempty(trialave); trialave = 1; end
    
    
    
    % Check if Chronux is loaded
    if ~exist('mtspectrumc','file');
        warning([' Chronux does not appear to be loaded. ' ...
            ' Make sure that Chronux' ...
            ' is installed and is in your Matlab path.']); end
    


    %params.tapers = get_tapers;
    %movingwind = [1 0.2];   
    Nshift = round(Nwind*(1-fract_overlap));

    Xstarts = 1:Nshift:(Nlfp-Nwind+1);
    Xstops = Xstarts + Nwind - 1;

    clear P
    
    % Get return dimensions
    i=1;
    Xsect = cellfunu(@(x) x(Xstarts(i):Xstops(i),:),X);
    [temp{1:nout-1}]=fname(Xsect{:});
    
    if isvector(temp); temp=temp(:); end
    sz = size(temp);
    
        % This assumes temp is a cell array. Might cause problems and need to
        % revert to earlier version; have not tested.
    OUTS = repmat({[]},1,nout);
    for j = 1:nout-1
        sz = size(temp{j});
        ending_dimension = length(sz);
        
        
        % Convert to GPU array if needed
        if is_gpuArray(X{1})
            OUTS{j} = zeros([sz,ones(1,10-ending_dimension),length(Xstarts)],'gpuArray');      % Store spectrogram T stuff in 11th dimension.
        else
            OUTS{j} = zeros([sz,ones(1,10-ending_dimension),length(Xstarts)]);      % Store spectrogram T stuff in 11th dimension.
        end
    end
    
    
    % Loop through all Xstarts
    for i = 1:length(Xstarts)
        Xsect = cellfunu(@(x) x(Xstarts(i):Xstops(i),:),X);
        [curr_outs{1:nout-1}]=fname(Xsect{:});           % Last arg out is for holding T, time.
        curr_outs = myvectorize(curr_outs);
        
        % Pack outputs
        for j = 1:nout-1
            %ending_dimension = length(size(curr_outs{j}));
            %OUTS{j} = cat(ending_dimension+1,OUTS{j},curr_outs{j});
            OUTS{j}(:,:,:,:,:,:,:,:,:,:,i) = curr_outs{j};
            
        end
        
        % Old cell array code for packing outputs
        %OUTS(1:nout-1) = cellfunu(@(X,Y) cat(3,X,Y),OUTS(1:nout-1),curr_outs);
    end
    
    for j = 1:nout-1;
        OUTS{j} = squeeze(OUTS{j});
    end

    %OUTS = cellfunu(@(X) permute(X,[1,3,2]),OUTS);
    if trialave; OUTS = cellfunu(@(X) mean(X,3),OUTS); end

    t = 1:Nlfp; t = t/fs;
    T = t(round(mean([Xstarts;Xstops])));
    T = T(:);
    
    OUTS{end} = T;
    
    if plot_on
        S=OUTS{1};
        F=OUTS{indexF}(:,1);
        
        if ndims(S) == 3
            figure
            for i = 1:size(S,3)
                clf; imagesc([min(T),max(T)],[min(F),max(F)],(S(:,:,i)));set(gca,'YDir','normal'); pause(1);
            end
        elseif ndims(S) < 3
            figure; imagesc([min(T),max(T)],[min(F),max(F)],(S));set(gca,'YDir','normal')
        else
            warning('Something wrong: too many dims of S');
        end
    end
    
    varargout = OUTS;
    clear OUTS;
end


function outs = myvectorize(outs)
    for i = 1:length(outs)
        if isvector(outs{i}); outs{i} = outs{i}(:); end
    end
end
