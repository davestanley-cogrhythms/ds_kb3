

function [P, f] = psd_wrapper(lfp,fs, mode)

    %if ~exist('N','var'); N = length(lfp); end
    if ~exist('fs','var'); fs = 1; end
    if ~exist('mode','var'); mode = 2; end          % Use Chronux mode by defualt
    
    
    if isvector(lfp)
        lfp=lfp(:);
    end
    
    N = size(lfp,1);
    Nwinds = size(lfp,2);
    
    switch mode
        case 1
            [Ptemp f] = pwelch(lfp(:,1),N,1,N,fs);
            ln = length(Ptemp); clear Ptemp     % Test run to get length
            P = zeros(ln,Nwinds);
            % Pwelch with only 1 window. We'll be averaging over trials anyways.
            for j = 1:Nwinds
                [P(:,j) f] = pwelch(lfp(:,j),N,1,N,fs); 
            end
            %P = mean(P,2);
    
        case 2
            params.Fs = fs;
            %params.trialave = 0;
            params.tapers = get_tapers;
            %tic
            [P f] = mtspectrumc(lfp,params); 
            %toc
    end

end