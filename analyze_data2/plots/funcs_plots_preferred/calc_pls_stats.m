
function [pls_stats, freqband_stats] = calc_pls_stats(f,pls,freqband_stats,varargin)

    p = inputParser;
    addOptional(p,'do_mean_ctgs',0,@isnumeric);     % If a 3rd dimension is present (ctgs) average across this; otherwise, leave this dimension present
    addOptional(p,'do_mean_freq',0,@isnumeric);     % Controls how freqbands are identified. If 0, take a single point at mean(freqband_stats); bandwidth determined by tapers. If 1, average all frequencies within freqband_stats (probably improper).
    parse(p,varargin{:});
    vars_pull(p.Results)
    
    if length(f) ~= size(pls,1)
        warning('f and pls dont match in size');
        return
    end

    if do_mean_freq
        index = f >= freqband_stats(1) & f <= freqband_stats(2);
    else
        %index = find(f >= mean(freqband_stats),1,'first');
        index = find_closest(f,mean(freqband_stats));
    end
    
    if do_mean_ctgs
        pls_stats = mean(mean(pls(index,:,:,:,:,:,:),3),1);
    else
        pls_stats = mean(pls(index,:,:,:,:,:,:),1);
        pls_stats = squeeze(pls_stats);
    end
end

