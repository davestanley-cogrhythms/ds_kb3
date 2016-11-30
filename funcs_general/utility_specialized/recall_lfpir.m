

function lfpir = recall_lfpir(lfp_sample_indices)
    
    % Get relative lfp indices
    
    lfpir = (lfp_sample_indices(1,1):lfp_sample_indices(1,3)) - lfp_sample_indices(1,2);

end