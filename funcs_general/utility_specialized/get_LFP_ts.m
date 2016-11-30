


function [currlfp,lfpia_stage] = get_LFP_ts(lfp_sample,curr_stage,md,lfp_sample_indices)

% [currlfp,lfpia_stage] = get_LFP_ts(lfp_sample,curr_stage,md,lfp_sample_indices)
% Gets LFP time series for each electrode. Does not do matching to the unit number.
% For this, should instead use get_unitLFP_ts
% 
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
%     lfpia_stage - Absolute lfp indices for the extracted 

    plot_on = 0;

    
    % % % % % % % % LFP % % % % % % % % % % 
    %unit_ind = md.unit_ind;
    
    % Get indices to pull out desired stage of data
    lfpir = recall_lfpir(lfp_sample_indices);
    boundsir = get_stagesir(curr_stage);
    lfpilr_stage = lfpir>=boundsir(1) & lfpir<=boundsir(2);         % Index logical relative (ilr)
    %lfpir_stage_old = recall_lfp_ind_stage(lfp_sample_indices,curr_stage);     % % Testing code
    
    currlfp = lfp_sample(lfpilr_stage,:,:);
    
    lfpia_stage = recall_lfpia(lfp_sample_indices);
    lfpia_stage = lfpia_stage(lfpilr_stage,:);
    
    if plot_on
        % Not really informative plotting...
        n=1:size(lfp_sample,1);
        plot(n(lfpilr_stage),currlfp);
        
    end
    
end