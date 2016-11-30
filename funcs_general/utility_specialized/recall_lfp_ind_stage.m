
% Function depreciated. Just use recall_lfpir and get_stagesir instead. See
% function get_lfp_and_spike_pairs.m for example code.
function index = recall_lfp_ind_stage(lfp_sample_indices,stage)

    % Pull out indices only for a specific stage of the sample_delay matrix
    % stage = 1 - pre sample on
    % stage = 2 - sample stage
    % stage = 3 - delay stage
    % stage = 4 - sample + delay stage
    
    dt = get_dt;
    
    lfp_indices = recall_lfpir(lfp_sample_indices);
    
    switch stage
        case 1
            index = lfp_indices < 0;
        case 2
            index = lfp_indices >= 0 & lfp_indices < round(1/dt);
        case 3
            index = lfp_indices >= round(1/dt);
        case 4
            index = lfp_indices >= 0;
    end
            
end