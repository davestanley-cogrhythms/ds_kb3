
function [pls_stats, freqband, timeband] = calc_pls3D_stats(f,f2,pls,freqband,timeband,varargin)


    p = inputParser;
    addOptional(p,'do_mean_ctgs',0,@isnumeric);     % If a 3rd dimension is present (ctgs) average across this; otherwise, leave this dimension present
    addOptional(p,'do_mean_freq',0,@isnumeric);     % Controls how freqbands are identified. If 0, take a single point at mean(freqband_stats); bandwidth determined by tapers. If 1, average all frequencies within freqband_stats (probably improper).
    addOptional(p,'ctg_singleton',0,@isnumeric);
    parse(p,varargin{:});
    
    do_mean_ctgs = p.Results.do_mean_ctgs;
    do_mean_freq = p.Results.do_mean_freq;
    ctg_singleton = p.Results.ctg_singleton;        % If there is only one ctg present and this has been squeezed out of the data.
                                                    % **Set to 1 when using group.data to get group.data stats.**
                                                    

    if ctg_singleton == 1
        time_dimension = 3;
    else
        time_dimension = 4;
    end
    
    % Check to make sure no empty data
    if isempty(f); error('f not assigned'); end
    if isempty(freqband); error('freqband not assigned'); end
    if isempty(f2); error('f2 not assigned'); end
    if isempty(timeband); error('timeband not assigned'); end
    
    % Check to make sure dimensions match
    if length(f) ~= size(pls,1)
        error('f and pls freq dimension dont match');
        return
    end
    if length(f2) ~= size(pls,time_dimension)
        error('f2 and pls time dimension dont match');
        return
    end
    
    % Calculate index and index2, along which we will slice the data.
    if do_mean_freq
        if length(freqband) < 2; error('To do mean, freqband must be lenght 2'); end
        if length(timeband) < 2; error('To do mean, timeband must be lenght 2'); end
        index = f >= freqband(1) & f <= freqband(2);
        index2 = f2 >= timeband(1) & f <= timeband(2);
    else
        index = find_closest(f,mean(freqband));
        index2 = find_closest(f2,mean(timeband));
    end
    
    % Slice data
    if ctg_singleton == 1
        pls_stats = mean(mean(pls(index,:,index2,:,:,:,:),1),3);    % Slice data along dims 1 and 3. Mean isn't necessary if not doing mean freqband
        pls_stats = squeeze(pls_stats);
    else
        pls_stats = mean(mean(pls(index,:,:,index2,:,:,:,:),1),4);    % Slice data along dims 1 and 4. Mean isn't necessary if not doing mean freqband
        if do_mean_ctgs
            pls_stats = mean(pls_stats,3);                           % Do mean along ctgs
        end
        pls_stats = squeeze(pls_stats);
    end
    
end

