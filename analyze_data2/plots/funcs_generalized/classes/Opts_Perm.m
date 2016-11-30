
classdef Opts_Perm
    % Class containing default settings
    properties
        % Perm switches    
        do_bh0 = 1;
        do_phi = 0;
        split_plusminus = 0;
        alpha0 = 0.05
        alpha_bh0 = 0.20;
        do_quantiles_mode = 0;
            chosen_quantile = 0.75;
            upper_quantile = 1;
        timeband_perm = [];
        tf_label_perm = 'DefaultPerm';
    end
end
    

