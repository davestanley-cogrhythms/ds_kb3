

function [lfp_all0] = load_lfp(fname,curr_stage,md)
    % function [coefs1,coefs2,data1,data2] = precalc_fft(fname,sfc_mode,curr_stage,filenum,i0,j0)
    
   
    % Load LFP data if needed
    load (fullfile(getpath('path_lfp_sample_delay'),fname));
    lfp_sample = double(lfp_sample);
    lfp_sample_indices = double(lfp_sample_indices); % Necessary for parfor, so knows data is loaded
    [lfp_all0,lfpia_stage] = get_LFP_ts(lfp_sample,curr_stage,md,lfp_sample_indices);

    
end
