

function group = group_normalize_specgram_to_baseline_time(group,specgram_baseline_time,normalize_within_elects, do_log)

    if ~exist('do_log','var')
        do_log = 0;
    end

    %normalize_within_elects = 1;   % 0 - normalize by mean of the pre-cue period.
                                    % 1 - normalize by each electrode individually. This within-electrode
                                    %    normalizing removes electrode-electrode variability.
                                    
    normalize_across_timerange = 1;
        normalizing_timerange = 0.1;

    for i = 1:length(group)
        if ~normalize_across_timerange
            ind = find(group(i).xdata2 >= specgram_baseline_time, 1, 'first');     % Find index to normalize along
        else
            ind = find( (group(i).xdata2 >= specgram_baseline_time - normalizing_timerange) & (group(i).xdata2 <= specgram_baseline_time + normalizing_timerange)); 
        end
        data = group(i).data;
        if ~normalize_within_elects
            data = data ./ repmat(mean(mean(data(:,:,ind),3),2),[1,size(data,2),size(data,3)]);
        else
            data = data ./ repmat(mean(data(:,:,ind),3),[1,1,size(data,3)]);
        end

        if do_log
            data = log(data);
        end
        
        group(i).data = data;
    
    end

end
