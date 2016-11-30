

function [avg_spt, spike_rate,Ntraces_below_thresh,fract_below_thresh] = format_spikes_per_trial(spikes_per_trial,curr_stage)


    %spt_threshold = get_spt_threshold; % Threshold for determining the minimum acceptable number of spikes per trial
    spt_threshold = 5;
    T = (diff(get_stagesir(curr_stage)) + 1)*get_dt; % Time duration of current stage (seconds)
    
    
    [avg_spt, spike_rate,Ntraces_below_thresh,fract_below_thresh] =  ... 
        cellfun(@(x) analyse_spt_stats(x,T,spt_threshold), spikes_per_trial);

end


function [avg_spt, spike_rate,Ntraces_below_thresh,fract_below_thresh] = format_spikes_per_trial(curr_spt,T,spt_threshold)


    Ntrials = length(curr_spt);
    tot_spikes = sum(curr_spt);
    avg_spt = sum(curr_spt) / Ntrials; % Average # spikes per trial
    spike_rate = avg_spt / T;

    below_thresh = curr_spt < spt_threshold;
    Ntraces_below_thresh = sum(below_thresh);
    fract_below_thresh = Ntraces_below_thresh / Ntrials;
    
end