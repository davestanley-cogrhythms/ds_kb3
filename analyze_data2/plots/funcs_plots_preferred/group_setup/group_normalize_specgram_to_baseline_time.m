

function group = group_normalize_specgram_to_baseline_time(group,specgram_baseline_time,normalize_within_elects)


    %normalize_within_elects = 1;   % 0 - normalize by mean of the pre-cue period.
                                    % 1 - normalize by each electrode individually. This within-electrode
                                    %    normalizing removes electrode-electrode variability.

    for i = 1:length(group)
        ind = find(group(i).xdata2 >= specgram_baseline_time, 1, 'first');     % Find index to normalize along
        data = group(i).data;
        if ~normalize_within_elects
            data = data ./ repmat(mean(data(:,:,ind),2),[1,size(data,2),size(data,3)]);
        else
            data = data ./ repmat(data(:,:,ind),[1,1,size(data,3)]);
        end

        data = log(data);
        
        group(i).data = data;
    
    end

end
