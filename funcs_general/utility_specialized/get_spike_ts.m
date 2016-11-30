


function [currspike] = get_spike_ts(curr_unit,curr_stage,md)

%     [currspike] = get_spike_ts(curr_unit,curr_stage,md)
%     INPUT
%     unit_ind - Units indices
%     curr_unit - Current unit to look at. The code will automagically pair this up with
%         the correct electrode based on metadata
%     curr_stage - Current stage of recording to look at (sample, delay, both, etc)
%     md - Metadata structure
% 
%     OUTPUT
%     currspike - Time series for spikes

    % % % % % % % % LFP % % % % % % % % % % 
    unit_ind = md.unit_ind;
        
    % % % % % % % % Spikes % % % % % % % % % % 
    boundsia = get_stagesia(curr_stage,md.sample_on_ind);
    currspike = get_spikes_n_range(unit_ind{curr_unit},[boundsia]');
    
end