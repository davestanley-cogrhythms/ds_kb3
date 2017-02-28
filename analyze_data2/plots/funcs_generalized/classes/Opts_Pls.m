
classdef Opts_Pls
    % Class containing default settings
    properties
        
        % Pls switches
        collapse_pls_to_days = 0;
        plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
        permdat2pls = 1;                % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
        perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
            perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
            perm2pls_dophi = 0;
            perm2pls_return_mode = 4;        % Return mode of perm2pls (4=zscore) (see perm2pls_swap)
            perm2pls_allow_signed = 0;       % If 1, alows perm2pls return to be +ve or -ve; else always positive. (see perm2pls_swap)
            perm2pls_split_plusminus = 0;
            sort_pls = 0;               % Sort into preferred and non-preferred
            swap_pls = [];
            do_diff = 1;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
            do_diff_percent = 0;
            do_abs_diff = 0;            % Take absolute value after doing diff.
            target_pls_format = 0;
        timeband_stats = [0.9 1.7];     % Delay by default.
        tf_label_stats = 'Default_Label';
        spectrogram2spectra_timeslice = 0;   % If working with a spectrogram, take a slice at time given by timeband_stats.
        spectrogram2ts_freqslice = 0;   % If working with a spectrogram, take a slice at time given by timeband_stats.
        spectrogram_normalize_to_baseline = 1;
            spectrogram_baseline_time = -0.75;       % During pre-cue
            spectrogram_baseline_dolog = 1;
        
    end
end
    

    

