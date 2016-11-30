
% Main function
function [spike_all0,lfp_all0] = preload_and_save(os,md)
    % Preload and store in variables.
    
    % Unpack os
    fname = os.fname;
    curr_stage = os.curr_stage;
    Nunits = length(md.unit_names);
    
    % Initialize
    spike_all0 = [];
    lfp_all0 = [];
    T_spk = [];
    T_lfp = [];
    T = [];
    
    
    % If not in spectrogram mode, just load the data normally.
    if ~os.is_spectrogram
        
        if os.needed_lfp
            lfp_all0 = load_lfp(fname,curr_stage,md);
        end
        
        if os.needed_spikes
            spike_all0 = get_spike_ts_all(Nunits,curr_stage,md);
            
            if strcmp(os.fname,'L112106.mat') && os.needed_lfp  % Shorten spikes to match length of LFP due to missing data
                spike_all0 = spike_all0(:,1:size(lfp_all0,2),:);
            end
        end
        
        
    end
    
end




