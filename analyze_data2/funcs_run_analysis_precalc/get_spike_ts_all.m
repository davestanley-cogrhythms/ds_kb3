
function currspike_all = get_spike_ts_all(Nunits,curr_stage,md)
    % function currspike_all = get_spike_ts_all(Nunits,curr_stage,md)

    j = 1;
    [currspike] = get_spike_ts(j,curr_stage,md);
    sz = size(currspike);
    currspike_all = zeros([sz, Nunits]);
    currspike_all(:,:,j) = currspike;
    
    for j = 2:Nunits
        currspike_all(:,:,j) = get_spike_ts(j,curr_stage,md);
    end
    
end
