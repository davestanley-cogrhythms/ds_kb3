

% FUNCTION DEPRECIATED. SEE THE FOLLOWING FUNCTIONS
% [currlfp lfpia_stage] = get_unitLFP_ts(lfp_sample,curr_unit,curr_stage,md,lfp_sample_indices);
% [currspike] = get_spike_ts(curr_unit,curr_stage,md);

function [currlfp currspike lfpia_stage] = get_lfp_and_spike_pairs(lfp_sample,curr_unit,curr_stage,md,lfp_sample_indices)

%     INPUT
%     lfp_sample - Full LFP matrix centered on sample on times
%     unit_ind - Units indices
%     curr_unit - Current unit to look at. The code will automagically pair this up with
%         the correct electrode based on metadata
%     curr_stage - Current stage of recording to look at (sample, delay, both, etc)
%     md - Metadata structure
%     lfp_sample_indices - LFP indicies. This should be packed alongside lfp_sample
% 
%     OUTPUT
%     currlfp - Time series for lfp
%     currspike - Time series for spikes
%     lfpia_stage - Absolute lfp indices for the extracted 

    plot_on = 0;

    
    % % % % % % % % LFP % % % % % % % % % % 
    curr_electrode = md.unit_LFP_pairs(2,curr_unit);
    unit_ind = md.unit_ind;
    
    % Get indices to pull out desired stage of data
    lfpir = recall_lfpir(lfp_sample_indices);
    boundsir = get_stagesir(curr_stage);
    lfpilr_stage = lfpir>=boundsir(1) & lfpir<=boundsir(2);         % Index logical relative (ilr)
    %lfpir_stage_old = recall_lfp_ind_stage(lfp_sample_indices,curr_stage);     % % Testing code
    
    currlfp = lfp_sample(lfpilr_stage,:,curr_electrode);
    
    lfpia_stage = recall_lfpia(lfp_sample_indices);
    lfpia_stage = lfpia_stage(lfpilr_stage,:);
    
    
    % % % % % % % % Spikes % % % % % % % % % % 
    
    % More testing code
%     boundsia2 = [lfp_sample_indices(:,2) lfp_sample_indices(:,3)];
%     sum(boundsia ~=boundsia2)
%     sia2 = get_stagesia(1,md.sample_on_ind);
%     sia3 = get_stagesia(3,md.sample_on_ind);
%     A = [sia2(:,1) md.sample_on_ind sia3(:,2)];
%     sum(lfp_sample_indices == A)
    
    boundsia = get_stagesia(curr_stage,md.sample_on_ind);
    currspike = get_spikes_n_range(unit_ind{curr_unit},[boundsia]');
    
    if plot_on
        % Not really informative plotting...
        n=1:size(lfp_sample,1);
        plot(n(lfpilr_stage),currlfp);
        
    end
    
end