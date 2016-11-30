

function [bad_lowfiring] = get_lowfiring(spikerate)

    % Collect spiking statistics to filter out neurons with spiking rates that
    % are too low to get reliable SFC curves
    %[avg_spt, spike_rate,Ntraces_below_thresh,fract_below_thresh] = format_spikes_per_trial(spikes_per_trial,curr_stage);
    %percent_below_thresh = N_below_thresh./Ntraces*100;
    
    %bad_lowfiring = percent_below_thresh(:,9) > 75; % Greater than 33% of trials have low firing rates! Proably a bad metric...
    bad_lowfiring = spikerate < 2.5;
    %bad_lowfiring = spikerate < 1;
    bad_lowfiring = bad_lowfiring(:)';

    
end

