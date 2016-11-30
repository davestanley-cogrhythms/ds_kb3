

function [numspikes] = get_numspikes_fromfile(curr_unit,curr_stage,md)
    % This actually returns the spike rate per second!

    unit_ind = md.unit_ind;
    boundsia = get_stagesia(curr_stage,md.sample_on_ind);
    currspike = get_spikes_n_range(unit_ind{curr_unit},[boundsia]');
    numspikes = mean(currspike);
    
    %duration = (mode(diff(boundsia'))+1)*get_dt;
    numspikes = numspikes / get_dt;       % Make it numspikes per second

end