

% Outdated functoin. Should use bounds_to_indices instead!

function sampleon_indices = recall_lfpia(lfp_sample_indices)

%     fprintf('Function outdated. \n Should remove. \n');
%     
%     sampleon_indices = (lfp_sample_indices(1,1):lfp_sample_indices(1,3)) - lfp_sample_indices(1,2);
%     sampleon_indices = repmat(sampleon_indices,size(lfp_sample_indices,1),1) + repmat(lfp_sample_indices(:,2),1,size(sampleon_indices,2));
%     sampleon_indices=sampleon_indices';
%     
    
    sampleon_indices = bounds_to_indices(lfp_sample_indices(:,[1 3])');

end