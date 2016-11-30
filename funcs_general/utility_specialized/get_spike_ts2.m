


function [spike_ts] = get_spike_ts2(unit_ind_cellarray,md)

%     INPUT
%     unit_ind - Units indices
%     curr_unit - Current unit to look at. The code will automagically pair this up with
%         the correct electrode based on metadata
%     curr_stage - Current stage of recording to look at (sample, delay, both, etc)
%     md - Metadata structure
% 
%     OUTPUT
%     currspike - Time series for spikes
    if iscell (unit_ind_cellarray)
        spike_ts = cellfun(@(x) spike_ind_to_ts_single (x,md),unit_ind_cellarray);
    elseif isa(unit_ind_cellarray,'double')
        spike_ts = spike_ind_to_ts_single (unit_ind_cellarray,md);
    end

    
        
end

function spike_ts = spike_ind_to_ts_single (unit_ind,md)

%     INPUT
%     unit_ind - Units indices
%     curr_unit - Current unit to look at. The code will automagically pair this up with
%         the correct electrode based on metadata
%     curr_stage - Current stage of recording to look at (sample, delay, both, etc)
%     md - Metadata structure
% 
%     OUTPUT
%     currspike - Time series for spikes
        
    % % % % % % % % Spikes % % % % % % % % % % 
    curr_stage = 5;     % Take all data
    boundsia = get_stagesia(curr_stage,md.sample_on_ind);
    spike_ts = get_spikes_n_range(unit_ind,[boundsia]');    
end


