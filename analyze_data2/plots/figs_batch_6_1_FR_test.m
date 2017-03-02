
%%
error('This is meant to run in cell mode. Dont F5!');


%% All figures

% % Clear
clearvars -except wrkspc_buffer

addpath('figs_batch_6_0_supporting');
addpath(genpath('funcs_plots_preferred'));
addpath(genpath('funcs_generalized'));

% % Load defaults
run_setdefaultfig

% % Create wrkspc_buffer structure if not already present
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end


% % Case 23; Ctgsetli Mode 0; - FFC
% g.mode_23_0 =  23.4013111011;
% g.mode_23_0 =  23.4013111013;
g.mode_23_0 =  23.4018111011;
% g.mode_23_0 =  23.4017101013;

% % Case 23; Ctgsetli Mode 40; - FFC
g.mode_23_40 =  23.4018103041;

% % Case 52; Ctgsetli Mode 0; - Units
g.mode_52_0 = 52.700001001;    
g.mode_52_2 = 52.700101001;     % Nwind = 300
% g.mode_52_2 = 52.700201001;     % Nwind = 600

% % Case 45; Ctgsetli Mode 0 - PSD
g.mode_45_0 = 45.6013111011;
% g.mode_45_0 = 45.6018111011;
% g.mode_45_0 = 45.6017101013;

% % Case 45; Ctgsetli Mode 40 - PSD
g.mode_45_40 = 45.6018103043;


%% Figure 6_10a - Plot unit FR based on bias - raw
clearvars -except wrkspc_buffer fv fv3 g
currfigname = 'Fg6_10a';

    % Setup params 
    clear group group0

    data_mode = 52.2;
    switch data_mode
        case 2.0          % SFC
            s.sfc_mode =  2.0511101;
            s.perm_mode = 52.700001001;
        case 2.2          % SFC
            %s.sfc_mode =  2.2511101;
            s.sfc_mode =  2.201711101;
            s.perm_mode = 52.700001001;
        case 2.3          % SFC all pairs
            s.sfc_mode =  2.301711101;
            s.perm_mode = 2.301711101;
        case 3          % SFC spectrogram
            s.sfc_mode =  3.201511101;
            s.perm_mode = 2.201511101;
        case 22         % FFC
            s.sfc_mode =  22.401811101;
            s.perm_mode = s.sfc_mode;
%             s.sfc_mode =  22.401511111;     % Partial FFC
%             s.perm_mode = s.sfc_mode;
%             s.sfc_mode =  22.451411103;     % FFC boundary vs non-boundary
            s.perm_mode = 52.700001001;
        case 23         % FFC spectrogram
            s.sfc_mode = g.mode_23_0;
            s.perm_mode = s.sfc_mode
        case 41         % PSD
            s.sfc_mode =  41.601811101;
            s.perm_mode = 41.601811101;
            %s.perm_mode = 52.700001001;         % Units
        case 45         % PSD spectrogram
            s.sfc_mode =  g.mode_45_0
            s.perm_mode = s.sfc_mode;
        case 52         % Units time series
            s.sfc_mode  = 52.700001001;
            s.perm_mode = g.mode_52_0
        case 52.2         % Units
            %s.sfc_mode  =  52.700201001;
            s.sfc_mode  =  g.mode_52_2;
            s.perm_mode =  g.mode_52_0;
            %s.sfc_mode =  41.601811101;
    end

    % Stage selection
    s.curr_stage_sfc = 3;
    s.curr_stage_sp = 3;
    if get_is_stage4(s.sfc_mode); s.curr_stage_sfc = 4; end
    if get_is_stage4(s.perm_mode); s.curr_stage_sp = 4; end


    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 0;
    opts_exclude.excludeL = 0;
    opts_exclude.excludeO = 1; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % Pls switches
    s.freqband_stats = [16 20];
    s.freqband_perm = [16 20];
%     s.freqband_stats = [10 12];
%     s.freqband_perm = [10 12];
%     s.freqband_stats = [20 30];
%     s.freqband_perm = [20 30];
    s.timeband_stats = [.6]; s.timeband_perm = [.6];
    s.tf_label_stats = 'Default'; s.tf_label_perm = 'Default';
    [tf_avail] = get_freqband_timeband(s.perm_mode,opts_exclude); s.tf_avail = tf_avail;
    
%     i=6; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
%     i=6; s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
    

    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 1;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % 0=Abs(Diff); 1=Diff;
        opts_pls.perm2pls_split_plusminus = 0;             % 0-Return everything; 1-return pos+negative in separate columns; 2-return only positive (pls.*Ramp(Diff)); 3-return only negative (pls.*Ramp(-Diff)); 4-return only significant cells (doesn't really have a use now)
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        opts_pls.do_diff_percent = 0;
        opts_pls.do_abs_diff = 0;            % Take absolute value after doing diff.
    opts_pls.target_pls_format = 0; % Convert pls to match this format!
    opts_pls.collapse_pls_to_days = 0;
    opts_pls.spectrogram2spectra_timeslice = 0;   % If working with a spectrogram, take a slice at time given by timeband_stats.
    opts_pls.spectrogram2ts_freqslice = 0;
    opts_pls.spectrogram_normalize_to_baseline = 0;          % Normalize spectrograms to pre-cue data to a value of 1.0
            opts_pls.spectrogram_baseline_time = -1.1;       % During pre-cue
            opts_pls.spectrogram_baseline_dolog = 1;
    
    

    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 2;      % 0-Either; 1-Both; 2-Positive; 3-Negative
    opts_perm.alpha0 = 0.05;
    opts_perm.alpha_bh0 = 0.2;
    opts_perm.alpha_bh0 = 0.05;
    opts_perm.do_quantiles_mode = 0;
        opts_perm.chosen_quantile = .15;
        opts_perm.upper_quantile = 0;
    
    % Map sp parameters
    s.sp_threshold = 10;
    
    % Group options
    s.do_group_collapse_pls2days = 0;
    s.do_group_normalize_specgram_to_baseline_time = 0;
        s.normalize_within_elects = 1;
        s.specgram_baseline_time = -1.1; s.gnsbt_do_log = 0;
    
    % Groupmode
    s.groupmode = 2.3;   % 0-Use default grouping (all pairs, enumerate over ctgs);
                       % 1:4-Select various subgroups
                       % 5:6-Separate into days
        s.examine_Sch_based_on_animal = 1;          % For animal L, do Cat/Dog; for O do Goc/Tad

    s.swap_mode = 0;

    s.group_do_merge = 0;
        s.groupmerge_operation = 0;

% % Plot switches
    s.plot_on_spect = 1;
    s.plot_on_scatter = 0;
    s.plot_on_bargraph = 0;
    s.plot_on_histogram = 0;
    s.plot_on_real_vs_imag = 0;
    s.plot_on_cells_vs_days = 0;
    s.plot_on_electrode_locations = 0;
    s.plot_on_electrode_distances = 0;
    
    s.doing_cat_vs_dog = (opts_pls.perm2pls && opts_pls.perm2pls_allow_signed == 1) || ...
        (opts_pls.permdat2pls == 1 && opts_pls.do_diff == 1 && opts_pls.do_abs_diff == 0) || ...
        s.group_do_merge == 1;
    
    s.do_custom_colourmap = s.doing_cat_vs_dog;
    s.do_custom_colourmap = false;


% % Plotting options
    paperfig_mode = 1;
    s.opts_PM3Dcs.paperfig_mode=paperfig_mode;
    s.opts_PM3Dcs.stats_mode = 0;
    s.opts_PM3Dcs.do_subplots = 0;
        s.opts_PM3Dcs.max_subplots_per_fig = 16;
    % % Spectrogram plotting options
    s.opts_PM3Dsp.paperfig_mode=paperfig_mode;
    s.opts_PM3Dsp.symmetric_axes = s.doing_cat_vs_dog || opts_pls.spectrogram_normalize_to_baseline || s.do_group_normalize_specgram_to_baseline_time;
    s.opts_PM3Dsp.symmetric_axes = [];
    s.opts_PM3Dsp.uniform_zaxis_allgroups = 1;           % Makes z-axis the same for all groups plotted
    s.opts_PM3Dsp.do_subplots = 1;
        s.opts_PM3Dsp.max_subplots_per_fig = 16;
    s.opts_PM3Dsp.show_range_stats = 1;
    s.opts_PM3Dsp.show_range_perm = 1;
        % Overlay Options - Transparency & Contours
        s.PM3Dsp_overlay_opts.do_transparency = 0;
        s.PM3Dsp_overlay_opts.do_contours = 0;
            s.overlay_raw_contours = 0;               % Overlays contours showing raw (non-diffed) FFC values.
            s.swap_in_groupdata_contours = 0;         % Overlays contours showing the same data being plotted in spectrogram (taken from group.data)
            s.swap_in_grouppairs_merge_pvals = 0;     % Overlays contours showing p values
        s.PM3Dsp_overlay_opts.contour_nv = [];
        s.PM3Dsp_overlay_opts.contour_linespec = {'k.'};
        % Stats
        s.PM3Dsp_stats_opts.stats_displaymode = 0;    % 0-no stats; 1 transparency; 2-contours; 3-both (overwrites default overlay settings above)
        s.PM3Dsp_stats_opts.statsfunc = [];
        s.PM3Dsp_stats_opts.stats_comparison = [0];
        s.PM3Dsp_stats_opts.contours_alphas = [0.01 0.001 0.0001];
        s.PM3Dsp_stats_opts.transparency_alpha = [0.01];
    
    s.opts_PSC.paperfig_mode = paperfig_mode;
    s.opts_PSC.remove_dependent = 0;
    s.opts_PSC.do_stats = 1;
    s.opts_PSC.do_diagtest = 0;
    s.opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    s.opts_PSC.opts_diagtest.threshold = 1.0; 
    s.opts_PDays.do_subplots = 1;
        s.opts_PDays.max_subplots_per_fig = 16;
        

% % Animal L
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 

% % Sample Stage
s.curr_stage_sp = 2;
opts_perm.split_plusminus = 2;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
opts_perm.split_plusminus = 3;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
%%
% % Delay Stage
s.curr_stage_sp = 3;
opts_perm.split_plusminus = 2;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
opts_perm.split_plusminus = 3;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);

%%

% % Animal O
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 

% % Sample Stage
s.curr_stage_sp = 2;
opts_perm.split_plusminus = 2;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
opts_perm.split_plusminus = 3;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
%
% % Delay Stage
s.curr_stage_sp = 3;
opts_perm.split_plusminus = 2;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
opts_perm.split_plusminus = 3;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);




clear out1 out2



%% Figure 6_10b - Plot unit FR based on bias - diff (for testing only; perhaps delete!)
clearvars -except wrkspc_buffer fv fv3 g
currfigname = 'Fg6_10b';

    % Setup params 
    clear group group0

    data_mode = 52.2;
    switch data_mode
        case 2.0          % SFC
            s.sfc_mode =  2.0511101;
            s.perm_mode = 52.700001001;
        case 2.2          % SFC
            %s.sfc_mode =  2.2511101;
            s.sfc_mode =  2.201711101;
            s.perm_mode = 52.700001001;
        case 2.3          % SFC all pairs
            s.sfc_mode =  2.301711101;
            s.perm_mode = 2.301711101;
        case 3          % SFC spectrogram
            s.sfc_mode =  3.201511101;
            s.perm_mode = 2.201511101;
        case 22         % FFC
            s.sfc_mode =  22.401811101;
            s.perm_mode = s.sfc_mode;
%             s.sfc_mode =  22.401511111;     % Partial FFC
%             s.perm_mode = s.sfc_mode;
%             s.sfc_mode =  22.451411103;     % FFC boundary vs non-boundary
            s.perm_mode = 52.700001001;
        case 23         % FFC spectrogram
            s.sfc_mode = g.mode_23_0;
            s.perm_mode = s.sfc_mode
        case 41         % PSD
            s.sfc_mode =  41.601811101;
            s.perm_mode = 41.601811101;
            %s.perm_mode = 52.700001001;         % Units
        case 45         % PSD spectrogram
            s.sfc_mode =  g.mode_45_0
            s.perm_mode = s.sfc_mode;
        case 52         % Units time series
            s.sfc_mode  = 52.700001001;
            s.perm_mode = g.mode_52_0
        case 52.2         % Units
            %s.sfc_mode  =  52.700201001;
            s.sfc_mode  =  g.mode_52_2;
            s.perm_mode =  g.mode_52_0;
            %s.sfc_mode =  41.601811101;
    end

    % Stage selection
    s.curr_stage_sfc = 3;
    s.curr_stage_sp = 3;
    if get_is_stage4(s.sfc_mode); s.curr_stage_sfc = 4; end
    if get_is_stage4(s.perm_mode); s.curr_stage_sp = 4; end


    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 0;
    opts_exclude.excludeL = 0;
    opts_exclude.excludeO = 1; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % Pls switches
    s.freqband_stats = [16 20];
    s.freqband_perm = [16 20];
%     s.freqband_stats = [10 12];
%     s.freqband_perm = [10 12];
%     s.freqband_stats = [20 30];
%     s.freqband_perm = [20 30];
    s.timeband_stats = [.6]; s.timeband_perm = [.6];
    s.tf_label_stats = 'Default'; s.tf_label_perm = 'Default';
    [tf_avail] = get_freqband_timeband(s.perm_mode,opts_exclude); s.tf_avail = tf_avail;
    
%     i=6; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
%     i=6; s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
    

    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 1;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 1;                % 0=Abs(Diff); 1=Diff;
        opts_pls.perm2pls_split_plusminus = 0;             % 0-Return everything; 1-return pos+negative in separate columns; 2-return only positive (pls.*Ramp(Diff)); 3-return only negative (pls.*Ramp(-Diff)); 4-return only significant cells (doesn't really have a use now)
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        opts_pls.do_diff_percent = 0;
        opts_pls.do_abs_diff = 0;            % Take absolute value after doing diff.
    opts_pls.target_pls_format = 0; % Convert pls to match this format!
    opts_pls.collapse_pls_to_days = 0;
    opts_pls.spectrogram2spectra_timeslice = 0;   % If working with a spectrogram, take a slice at time given by timeband_stats.
    opts_pls.spectrogram2ts_freqslice = 0;
    opts_pls.spectrogram_normalize_to_baseline = 0;          % Normalize spectrograms to pre-cue data to a value of 1.0
            opts_pls.spectrogram_baseline_time = -1.1;       % During pre-cue
            opts_pls.spectrogram_baseline_dolog = 1;
    
    

    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;      % 0-Either; 1-Both; 2-Positive; 3-Negative
    opts_perm.alpha0 = 0.05;
    opts_perm.alpha_bh0 = 0.2;
    opts_perm.alpha_bh0 = 0.05;
    opts_perm.do_quantiles_mode = 0;
        opts_perm.chosen_quantile = .15;
        opts_perm.upper_quantile = 0;
    
    % Map sp parameters
    s.sp_threshold = 10;
    
    % Group options
    s.do_group_collapse_pls2days = 0;
    s.do_group_normalize_specgram_to_baseline_time = 0;
        s.normalize_within_elects = 1;
        s.specgram_baseline_time = -1.1; s.gnsbt_do_log = 0;
    
    % Groupmode
    s.groupmode = 2;   % 0-Use default grouping (all pairs, enumerate over ctgs);
                       % 1:4-Select various subgroups
                       % 5:6-Separate into days
        s.examine_Sch_based_on_animal = 1;          % For animal L, do Cat/Dog; for O do Goc/Tad

    s.swap_mode = 0;

    s.group_do_merge = 0;
        s.groupmerge_operation = 0;

% % Plot switches
    s.plot_on_spect = 1;
    s.plot_on_scatter = 0;
    s.plot_on_bargraph = 0;
    s.plot_on_histogram = 0;
    s.plot_on_real_vs_imag = 0;
    s.plot_on_cells_vs_days = 0;
    s.plot_on_electrode_locations = 0;
    s.plot_on_electrode_distances = 0;
    
    s.doing_cat_vs_dog = (opts_pls.perm2pls && opts_pls.perm2pls_allow_signed == 1) || ...
        (opts_pls.permdat2pls == 1 && opts_pls.do_diff == 1 && opts_pls.do_abs_diff == 0) || ...
        s.group_do_merge == 1;
    
    s.do_custom_colourmap = s.doing_cat_vs_dog;
    s.do_custom_colourmap = false;


% % Plotting options
    paperfig_mode = 1;
    s.opts_PM3Dcs.paperfig_mode=paperfig_mode;
    s.opts_PM3Dcs.stats_mode = 0;
    s.opts_PM3Dcs.do_subplots = 0;
        s.opts_PM3Dcs.max_subplots_per_fig = 16;
    % % Spectrogram plotting options
    s.opts_PM3Dsp.paperfig_mode=paperfig_mode;
    s.opts_PM3Dsp.symmetric_axes = s.doing_cat_vs_dog || opts_pls.spectrogram_normalize_to_baseline || s.do_group_normalize_specgram_to_baseline_time;
    s.opts_PM3Dsp.symmetric_axes = [];
    s.opts_PM3Dsp.uniform_zaxis_allgroups = 1;           % Makes z-axis the same for all groups plotted
    s.opts_PM3Dsp.do_subplots = 1;
        s.opts_PM3Dsp.max_subplots_per_fig = 16;
    s.opts_PM3Dsp.show_range_stats = 1;
    s.opts_PM3Dsp.show_range_perm = 1;
        % Overlay Options - Transparency & Contours
        s.PM3Dsp_overlay_opts.do_transparency = 0;
        s.PM3Dsp_overlay_opts.do_contours = 0;
            s.overlay_raw_contours = 0;               % Overlays contours showing raw (non-diffed) FFC values.
            s.swap_in_groupdata_contours = 0;         % Overlays contours showing the same data being plotted in spectrogram (taken from group.data)
            s.swap_in_grouppairs_merge_pvals = 0;     % Overlays contours showing p values
        s.PM3Dsp_overlay_opts.contour_nv = [];
        s.PM3Dsp_overlay_opts.contour_linespec = {'k.'};
        % Stats
        s.PM3Dsp_stats_opts.stats_displaymode = 0;    % 0-no stats; 1 transparency; 2-contours; 3-both (overwrites default overlay settings above)
        s.PM3Dsp_stats_opts.statsfunc = [];
        s.PM3Dsp_stats_opts.stats_comparison = [0];
        s.PM3Dsp_stats_opts.contours_alphas = [0.01 0.001 0.0001];
        s.PM3Dsp_stats_opts.transparency_alpha = [0.01];
    
    s.opts_PSC.paperfig_mode = paperfig_mode;
    s.opts_PSC.remove_dependent = 0;
    s.opts_PSC.do_stats = 1;
    s.opts_PSC.do_diagtest = 0;
    s.opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    s.opts_PSC.opts_diagtest.threshold = 1.0; 
    s.opts_PDays.do_subplots = 1;
        s.opts_PDays.max_subplots_per_fig = 16;
        

% % Animal L
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 

% % Sample Stage
s.curr_stage_sp = 2;
opts_perm.split_plusminus = 2;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
opts_perm.split_plusminus = 3;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
%%
% % Delay Stage
s.curr_stage_sp = 3;
opts_perm.split_plusminus = 2;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
opts_perm.split_plusminus = 3;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);



% % Animal O
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 

% % Sample Stage
s.curr_stage_sp = 2;
opts_perm.split_plusminus = 2;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
opts_perm.split_plusminus = 3;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);

% % Delay Stage
s.curr_stage_sp = 3;
opts_perm.split_plusminus = 2;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
opts_perm.split_plusminus = 3;      % 0-Either; 1-Both; 2-Positive; 3-Negative        
[wrkspc_buffer, out1] = Fg_6_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);




clear out1 out2


