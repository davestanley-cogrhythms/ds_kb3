

function [P, f] = psd_wrapper(varargin)
%     [P, f] = psd_wrapper(X,fs,'mode',mode,'params',params)
%     Calculates power spectral density using one of several
%     different methods, thus serving as a wrapper. 
%     FORMS
%         [P, f] = psd_wrapper(X)
%         [P, f] = psd_wrapper(X,options)
%     INPUTS
%         t - times (vector)
%         X - data (vector/matrix)
%         options - specifies options in the form of name and value pairs
%             (i.e. 'fs',1024)
%     OPTIONS
%         fs - sampling rate (Default=1)
%         mode - 1-Default Matlab PSD
%                2-Chronux (Default)
%         params - paramater structure to be fed into Chronux commands
%                       - note - in all cases, fs will override params.Fs.
%                       - only used if mode=2
%     OUTPUTS
%         [P, f] - power and frequency
% 
%     CONTACT: David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
% 


    % Pull inputs to local workspace
    p = inputParser;
    addRequired(p,'X',@isnumeric);
    addOptional(p,'fs',[],@isnumeric);
    addOptional(p,'mode',[],@isnumeric);     % Default mode is chronux
    addOptional(p,'params',[]);    
    parse(p,varargin{:});
    vars_pull(p.Results);
    mymode=p.Results.mode;
    
    % Set defaults
    if isempty(fs); fs = 1; end
    if isempty(mymode); mymode = 2; end
    
    % Setup X
    if isvector(X); X=X(:); end
    
    % Setup parameters
    params.Fs = fs;
    

%     %if ~exist('N','var'); N = length(lfp); end
%     if ~exist('mode','var'); mode = 2; end          % Use Chronux mode by defualt
    
    
    if isvector(X)
        X=X(:);
    end
    
    N = size(X,1);
    Nwinds = size(X,2);
    
    
    % Check if Chronux is loaded
    if mymode == 2 && ~exist('mtspectrumc','file');
        fprintf([' Error. Chronux does not appear to be loaded. Try switching to mode = 1 to use \n' ...
            ' default Matlab commands. Alternatively, make sure that Chronux \n' ...
            ' is installed and is in your Matlab path. Exiting. \n']); return; end
    
    
    switch mymode
        case 1
            [Ptemp f] = pwelch(X(:,1),N,1,N,fs);
            ln = length(Ptemp); clear Ptemp     % Test run to get length
            P = zeros(ln,Nwinds);
            % Pwelch with only 1 window. We'll be averaging over trials anyways.
            for j = 1:Nwinds
                [P(:,j), f] = pwelch(X(:,j),N,1,N,fs); 
            end
            %P = mean(P,2);
    
        case 2
            
            %tic
            [P, f] = mtspectrumc(X,params); 
            %toc
    end

end



function vars_pull(s)
    for n = fieldnames(s)'
        name = n{1};
        value = s.(name);
        assignin('caller',name,value);
    end
end