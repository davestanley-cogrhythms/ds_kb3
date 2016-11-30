
%%
error('This is meant to run in cell mode. Dont F5!');


%% All figures

% Clear
clearvars -except wrkspc_buffer fv fv3

addpath('figs_batch_5_0_supporting');
addpath(genpath('funcs_plots_preferred'));
addpath(genpath('funcs_generalized'));

% Load defaults
fv3 = figs_batch_defaults_3; 
run_setdefaultfig

% Create wrkspc_buffer structure if not already present
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end


%% Figure 5_00a - Analysis of average responses - spectrograms
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg5_00a';

% Setup params 
data_mode = 45;
switch data_mode
    case 2.2          % SFC
        s.sfc_mode =  2.201511101;
        s.perm_mode = 2.201511101;
    case 2.3          % SFC all pairs
        s.sfc_mode =  2.301711101;
        s.perm_mode = 2.301711101;
    case 3          % SFC spectrogram
        s.sfc_mode =  3.201511101;
        s.perm_mode = 2.201511101;
    case 22         % FFC
        s.sfc_mode =  22.401811101;
        s.perm_mode = 22.401811101;
    case 23         % FFC spectrogram
        s.sfc_mode =  23.401411101;
        s.perm_mode = 22.401811101;
    case 41         % PSD
        s.sfc_mode =  41.601811101;
        s.perm_mode = 41.601811101;
    case 45         % PSD spectrogram
        s.sfc_mode =  45.601311101;
        s.perm_mode = 41.601811101;
    case 52         % Units time series
        s.sfc_mode  = 52.700001001;
        s.perm_mode = 52.700001001;
    case 52.2         % Units
        s.sfc_mode  = 52.700201001;
        s.perm_mode = 52.700001001;
end

s.is_spectrogram = get_is_spectrogram(s.sfc_mode);

% Stage selection
s.curr_stage_sfc = 3;
s.curr_stage_sp = 3;

if get_is_stage4(s.sfc_mode); s.curr_stage_sfc = 4; end
if get_is_stage4(s.perm_mode); s.curr_stage_sp = 4; end

    % Pls switches
    s.freqband_stats = [16 20];
    s.freqband_stats_perm = [16 20];
%     s.freqband_stats = [10 12];
%     s.freqband_stats_perm = [10 12];

    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 0;
    opts_exclude.excludeL = 0;
    opts_exclude.excludeO = 1; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC

    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
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
        

    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;      % 0-Either; 1-Both; 2-Positive; 3-Negative
    opts_perm.alpha0 = 0.05;
    opts_perm.alpha_bh0 = 0.2;
%     opts_perm.alpha_bh0 = 0.05;
    
    % Map sp parameters
    s.sp_threshold = 10;
    
    % Group options
    s.do_group_collapse_pls2days = 0;
    
    
% Groupmode stuff    
s.groupmode = 0;   % 0-Use default grouping (all pairs, enumerate over ctgs);
                 % 1-Group based on permutation test output for diff.
                 % 2-Group based on permutation test output all.
    s.swap_mode = 0;

s.group_do_merge = 0;
    s.grouppmerge_do_percent = 1;

% % % % % % % % % % % % % % % % % % % % N Cells % % % % % % % % % % % % %
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 1;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % 0=Abs(Diff); 1=Diff;
        opts_pls.perm2pls_split_plusminus = 0;             % 0-Return everything; 1-return pos+negative in separate columns; 2-return only positive (pls.*Ramp(Diff)); 3-return only negative (pls.*Ramp(-Diff)); 4-return only significant cells (doesn't really have a use now)
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        opts_pls.do_diff_percent = 0;
        opts_pls.do_abs_diff = 0;            % Take absolute value after doing diff.

opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1;
opts_pls.perm2pls_split_plusminus = 0;
[wrkspc_buffer, out] = Fg_5_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)

opts_pls.perm2pls_split_plusminus = 2;
[wrkspc_buffer, out] = Fg_5_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)

opts_pls.perm2pls_split_plusminus = 3;
[wrkspc_buffer, out] = Fg_5_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)

% pause

opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0;
opts_pls.perm2pls_split_plusminus = 0;
[wrkspc_buffer, out] = Fg_5_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)

opts_pls.perm2pls_split_plusminus = 2;
[wrkspc_buffer, out] = Fg_5_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)

opts_pls.perm2pls_split_plusminus = 3;
[wrkspc_buffer, out] = Fg_5_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)

% pause

% % % % % % % % % % % % % % % % % % % % Z-scores Cells % % % % % % % % % % % % %
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
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

opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1;
opts_pls.perm2pls_allow_signed = 0;
[wrkspc_buffer, out] = Fg_5_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)

opts_pls.perm2pls_allow_signed = 1;
[wrkspc_buffer, out] = Fg_5_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)


% pause

opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0;
opts_pls.perm2pls_allow_signed = 0;
[wrkspc_buffer, out] = Fg_5_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)

opts_pls.perm2pls_allow_signed = 1;
[wrkspc_buffer, out] = Fg_5_00_generalized(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)



%% Figure 5_00b - Analysis of significant vs non-significant electrodes
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg5_00b';

% Setup params 
    clear group group0

    data_mode = 23;
    switch data_mode
        case 2.2          % SFC
            s.sfc_mode =  2.201511101;
            s.perm_mode = 2.201511101;
        case 2.3          % SFC all pairs
            s.sfc_mode =  2.301711101;
            s.perm_mode = 2.301711101;
        case 3          % SFC spectrogram
            s.sfc_mode =  3.201511101;
            s.perm_mode = 2.201511101;
        case 22         % FFC
            s.sfc_mode =  22.401811101;
            s.perm_mode = 22.401811101;
        case 23         % FFC spectrogram
            s.sfc_mode =  23.4013111011;
            s.perm_mode = s.sfc_mode;
        case 41         % PSD
            s.sfc_mode =  41.601811101;
            s.perm_mode = 41.601811101;
        case 45         % PSD spectrogram
            s.sfc_mode =  45.6013111011;
            s.perm_mode = s.sfc_mode;
        case 52         % Units time series
            s.sfc_mode  = 52.700001001;
            s.perm_mode = 52.700001001;
        case 52.2         % Units
            s.sfc_mode  = 52.700201000;
            s.perm_mode = 52.700001001;
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
    opts_exclude.excludeL = 1;
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % Pls switches
    s.freqband_stats = [16 20];
    s.freqband_perm = [16 20];
    s.timeband_stats = [.6]; s.timeband_perm = [.6];
    s.tf_label_stats = 'Default'; s.tf_label_perm = 'Default';
    [tf_avail] = get_freqband_timeband(s.sfc_mode,opts_exclude); s.tf_avail = tf_avail;
    
    i=1; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
    i=1; s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
    

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

    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;      % 0-Either; 1-Both; 2-Positive; 3-Negative
    opts_perm.alpha0 = 0.05;
    opts_perm.alpha_bh0 = 0.2;
    %opts_perm.alpha_bh0 = 0.05;
    opts_perm.do_quantiles_mode = 0;
        opts_perm.chosen_quantile = .75;
        opts_perm.upper_quantile = 1;
    
    % Map sp parameters
    s.sp_threshold = 10;
    
    % Group options
    s.do_group_collapse_pls2days = 0;
    
    % Groupmode
    s.groupmode = 0;   % 0-Use default grouping (all pairs, enumerate over ctgs);
                 % 1-Group based on permutation test output for diff.
                 % 2-Group based on permutation test output all.
    s.swap_mode = 0;

    s.group_do_merge = 0;
        s.grouppmerge_do_percent = 1;

% % Plot switches
    s.plot_on_spect = 1;
    s.plot_on_scatter = 0;
    s.plot_on_bargraph = 0;
    s.plot_on_histogram = 0;
    s.plot_on_real_vs_imag = 0;
    s.plot_on_cells_vs_days = 0;
    
    s.doing_cat_vs_dog = (opts_pls.perm2pls && opts_pls.perm2pls_allow_signed == 1) || ...
        (opts_pls.permdat2pls == 1 && opts_pls.do_diff == 1 && opts_pls.do_abs_diff == 0);
    
    s.do_custom_colourmap = s.doing_cat_vs_dog;
    s.do_custom_colourmap = false;


% % Plotting options
    paperfig_mode = 1;
    s.opts_PM3Dcs.paperfig_mode=paperfig_mode;
    s.opts_PM3Dcs.stats_mode = 0;
    s.opts_PM3Dcs.do_subplots = 0;
        s.opts_PM3Dcs.max_subplots_per_fig = 16;
    s.opts_PM3Dsp.paperfig_mode=paperfig_mode;
    s.opts_PM3Dsp.symmetric_axes = s.doing_cat_vs_dog;
    s.opts_PM3Dsp.uniform_zaxis_allgroups = 1;           % Makes z-axis the same for all groups plotted
    s.opts_PM3Dsp.do_subplots = 1;
        s.opts_PM3Dsp.max_subplots_per_fig = 16;
    s.opts_PSC.paperfig_mode = paperfig_mode;
    s.opts_PSC.remove_dependent = 0;
    s.opts_PSC.do_stats = 1;
    s.opts_PSC.do_diagtest = 0;
    s.opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    s.opts_PSC.opts_diagtest.threshold = 1.0; 
    s.opts_PDays.do_subplots = 1;
        s.opts_PDays.max_subplots_per_fig = 16;
        
%
        
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1 ; 
[tf_avail] = get_freqband_timeband(s.sfc_mode,opts_exclude); s.tf_avail = tf_avail;



i=3; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
i=6; s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
    

s.groupmode = 1
opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
opts_pls.spectrogram2spectra_timeslice = 0;
s.plot_on_spect = 1;
s.plot_on_bargraph = 0;
[wrkspc_buffer, out] = Fg_5_00b_generalized_updated(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)


s.groupmode = 1
opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
opts_pls.spectrogram2spectra_timeslice = 1;
s.plot_on_spect = 0;
s.plot_on_bargraph = 1;
[wrkspc_buffer, out] = Fg_5_00b_generalized_updated(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)

s.groupmode = 3
opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
opts_pls.spectrogram2spectra_timeslice = 1;
s.plot_on_spect = 1;
s.plot_on_bargraph = 0;
[wrkspc_buffer, out] = Fg_5_00b_generalized_updated(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)

%% Figure 5_00c - Analysis of positive vs negative responses
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg5_00c';

% Setup params 
    clear group group0

    data_mode = 23;
    switch data_mode
        case 2.2          % SFC
            s.sfc_mode =  2.201511101;
            s.perm_mode = 2.201511101;
        case 2.3          % SFC all pairs
            s.sfc_mode =  2.301711101;
            s.perm_mode = 2.301711101;
        case 3          % SFC spectrogram
            s.sfc_mode =  3.201511101;
            s.perm_mode = 2.201511101;
        case 22         % FFC
            s.sfc_mode =  22.401811101;
            s.perm_mode = 22.401811101;
        case 23         % FFC spectrogram
            s.sfc_mode =  23.4013111011;
            s.perm_mode = s.sfc_mode;
        case 41         % PSD
            s.sfc_mode =  41.601811101;
            s.perm_mode = 41.601811101;
        case 45         % PSD spectrogram
            s.sfc_mode =  45.6013111011;
            s.perm_mode = s.sfc_mode;
        case 52         % Units time series
            s.sfc_mode  = 52.700001001;
            s.perm_mode = 52.700001001;
        case 52.2         % Units
            s.sfc_mode  = 52.700201000;
            s.perm_mode = 52.700001001;
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
    opts_exclude.excludeL = 1;
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % Pls switches
    s.freqband_stats = [16 20];
    s.freqband_perm = [16 20];
    s.timeband_stats = [.6]; s.timeband_perm = [.6];
    s.tf_label_stats = 'Default'; s.tf_label_perm = 'Default';
    [tf_avail] = get_freqband_timeband(s.sfc_mode,opts_exclude); s.tf_avail = tf_avail;
    
    i=1; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
    i=1; s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
    

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

    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;      % 0-Either; 1-Both; 2-Positive; 3-Negative
    opts_perm.alpha0 = 0.05;
    opts_perm.alpha_bh0 = 0.2;
%     opts_perm.alpha_bh0 = 0.05;
    opts_perm.do_quantiles_mode = 0;
        opts_perm.chosen_quantile = .75;
        opts_perm.upper_quantile = 1;
    
    % Map sp parameters
    s.sp_threshold = 10;
    
    % Group options
    s.do_group_collapse_pls2days = 0;
    
    % Groupmode
    s.groupmode = 0;   % 0-Use default grouping (all pairs, enumerate over ctgs);
                 % 1-Group based on permutation test output for diff.
                 % 2-Group based on permutation test output all.
    s.swap_mode = 0;

    s.group_do_merge = 0;
        s.grouppmerge_do_percent = 1;

% % Plot switches
    s.plot_on_spect = 1;
    s.plot_on_scatter = 0;
    s.plot_on_bargraph = 0;
    s.plot_on_histogram = 0;
    s.plot_on_real_vs_imag = 0;
    s.plot_on_cells_vs_days = 0;
    
    s.doing_cat_vs_dog = (opts_pls.perm2pls && opts_pls.perm2pls_allow_signed == 1) || ...
        (opts_pls.permdat2pls == 1 && opts_pls.do_diff == 1 && opts_pls.do_abs_diff == 0);
    
    s.do_custom_colourmap = s.doing_cat_vs_dog;
    s.do_custom_colourmap = false;


% % Plotting options
    paperfig_mode = 1;
    s.opts_PM3Dcs.paperfig_mode=paperfig_mode;
    s.opts_PM3Dcs.stats_mode = 0;
    s.opts_PM3Dcs.do_subplots = 0;
        s.opts_PM3Dcs.max_subplots_per_fig = 16;
    s.opts_PM3Dsp.paperfig_mode=paperfig_mode;
    s.opts_PM3Dsp.symmetric_axes = s.doing_cat_vs_dog;
    s.opts_PM3Dsp.uniform_zaxis_allgroups = 1;           % Makes z-axis the same for all groups plotted
    s.opts_PM3Dsp.do_subplots = 1;
        s.opts_PM3Dsp.max_subplots_per_fig = 16;
    s.opts_PSC.paperfig_mode = paperfig_mode;
    s.opts_PSC.remove_dependent = 0;
    s.opts_PSC.do_stats = 1;
    s.opts_PSC.do_diagtest = 0;
    s.opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    s.opts_PSC.opts_diagtest.threshold = 1.0; 
    s.opts_PDays.do_subplots = 1;
        s.opts_PDays.max_subplots_per_fig = 16;
        
%
        
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 
[tf_avail] = get_freqband_timeband(s.sfc_mode,opts_exclude); s.tf_avail = tf_avail;



i=2; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
i=3; s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
    

s.groupmode = 1
opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
opts_pls.spectrogram2spectra_timeslice = 0;
s.plot_on_spect = 1;
s.plot_on_bargraph = 0;
opts_perm.split_plusminus = 2;      % 0-Either; 1-Both; 2-Positive; 3-Negative
[wrkspc_buffer, out] = Fg_5_00b_generalized_updated(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)


s.groupmode = 1
opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
opts_pls.spectrogram2spectra_timeslice = 0;
s.plot_on_spect = 1;
s.plot_on_bargraph = 0;
opts_perm.split_plusminus = 3;      % 0-Either; 1-Both; 2-Positive; 3-Negative
[wrkspc_buffer, out] = Fg_5_00b_generalized_updated(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)



%% Figure 5_01a - Analysis of unit biases grouped by FFC bias
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg5_01a';

% Setup params 
    clear group group0 s

    data_mode = 52.2;
    
            s.sfc_mode  =  52.700201001;
            s.perm_mode =  52.700001001;
            %s.sfc_mode =  41.601811101;
            s.perm_mode =  23.4013111011;
            s.perm_mode =  23.4018111013;
    

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
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % Pls switches
    s.freqband_stats = [16 20];
    s.freqband_perm = [16 20];
%     s.freqband_stats = [10 12];
%     s.freqband_perm = [10 12];
    s.timeband_stats = [.6]; s.timeband_perm = [.6];
    s.tf_label_stats = 'Default'; s.tf_label_perm = 'Default';
    [tf_avail] = get_freqband_timeband(s.perm_mode,opts_exclude); s.tf_avail = tf_avail;
    

    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 1;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 1;                % 0=Abs(Diff); 1=Diff;
        opts_pls.perm2pls_split_plusminus = 0;             % 0-Return everything; 1-return pos+negative in separate columns; 2-return only positive (pls.*Ramp(Diff)); 3-return only negative (pls.*Ramp(-Diff)); 4-return only significant cells (doesn't really have a use now)
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 1;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        opts_pls.do_diff_percent = 0;
        opts_pls.do_abs_diff = 0;            % Take absolute value after doing diff.
    opts_pls.target_pls_format = 0; % Convert pls to match this format!
    opts_pls.collapse_pls_to_days = 0;
    opts_pls.spectrogram2spectra_timeslice = 0;   % If working with a spectrogram, take a slice at time given by timeband_stats.
    opts_pls.spectrogram2ts_freqslice = 0;

    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;      % 0-Either; 1-Both; 2-Positive; 3-Negative
    opts_perm.alpha0 = 0.05;
    opts_perm.alpha_bh0 = 0.2;
    %opts_perm.alpha_bh0 = 0.05;
    opts_perm.do_quantiles_mode = 0;
        opts_perm.chosen_quantile = .15;
        opts_perm.upper_quantile = 0;
    
    % Map sp parameters
    s.sp_threshold = 40;
    
    % Group options
    s.do_group_collapse_pls2days = 0;
    
    % Groupmode
    s.groupmode = 1;   % 0-Use default grouping (all pairs, enumerate over ctgs);
                 % 1-Group based on permutation test output for diff.
                 % 2-Group based on permutation test output all.
    s.swap_mode = 0;

    s.group_do_merge = 0;
        s.grouppmerge_do_percent = 1;

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
        (opts_pls.permdat2pls == 1 && opts_pls.do_diff == 1 && opts_pls.do_abs_diff == 0);
    
    s.do_custom_colourmap = s.doing_cat_vs_dog;
    s.do_custom_colourmap = false;


% % Plotting options
    paperfig_mode = 1;
    s.opts_PM3Dcs.paperfig_mode=paperfig_mode;
    s.opts_PM3Dcs.stats_mode = 0;
    s.opts_PM3Dcs.do_subplots = 0;
        s.opts_PM3Dcs.max_subplots_per_fig = 16;
    s.opts_PM3Dsp.paperfig_mode=paperfig_mode;
    s.opts_PM3Dsp.symmetric_axes = s.doing_cat_vs_dog;
    s.opts_PM3Dsp.uniform_zaxis_allgroups = 1;           % Makes z-axis the same for all groups plotted
    s.opts_PM3Dsp.do_subplots = 1;
        s.opts_PM3Dsp.max_subplots_per_fig = 16;
    s.opts_PM3Dsp.show_range_stats = 1;
    s.opts_PM3Dsp.show_range_perm = 1;
    s.opts_PSC.paperfig_mode = paperfig_mode;
    s.opts_PSC.remove_dependent = 0;
    s.opts_PSC.do_stats = 1;
    s.opts_PSC.do_diagtest = 0;
    s.opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    s.opts_PSC.opts_diagtest.threshold = 1.0; 
    s.opts_PDays.do_subplots = 1;
        s.opts_PDays.max_subplots_per_fig = 16;
        
    do_animalL = 1;
    do_animalO = 1;
    
    do_positives = 1;
    do_negatives = 1;
    
    
    % Animal L
    if do_animalL
        opts_exclude.excludeL = 0;
        opts_exclude.excludeO = 1; 
        [tf_avail] = get_freqband_timeband(s.perm_mode,opts_exclude); s.tf_avail = tf_avail;

        % Positive bias (towards Cat/Goc)
        if do_positives
        %     i=1; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
        %          s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
        %     [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
        %     
        %     i=2; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
        %          s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
        %     [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);

            i=4; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                 s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);

%             i=5; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
%                  s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
%             [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
% 
%             i=6; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
%                  s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
%             [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
% 
        %     i=7; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
        %          s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
        %     [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
        end

        % Negative bias (towards Dog/Tad)
        if do_negatives
            i=3; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                 s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
        end
    end
        
    % Animal O, 
    if do_animalO
        opts_exclude.excludeL = 1;
        opts_exclude.excludeO = 0; 
        [tf_avail] = get_freqband_timeband(s.perm_mode,opts_exclude); s.tf_avail = tf_avail;

        % Positive bias (towards Cat/Goc; modes 1, 2, 8, 11)
        if do_positives
            %     i=2; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            %          s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            %     [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
            %     
            %     i=3; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            %          s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            %     [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
            %     
            %     i=8; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            %          s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            %     [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);

                i=11; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                      s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
        end

        % Negative bias (towards Dog/Tad; modes 3, 9, 10)
        if do_negatives
        %     i=1; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
        %          s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
        %     [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
% 
%             i=9; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
%                  s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
%             [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);

            i=10; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                 s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
        end
    end
    


%% Figure 5_01b - Analysis of unit bias grouped by POWER bias
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg5_01b';

% Setup params 
    clear group group0 s

    data_mode = 52.2;
    
            s.sfc_mode  =  52.700201001;
            s.perm_mode =  52.700001001;
            %s.sfc_mode =  41.601811101;
            s.perm_mode =  23.4013111011;
            s.perm_mode =  45.6013111010;
    

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
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % Pls switches
    s.freqband_stats = [16 20];
    s.freqband_perm = [16 20];
%     s.freqband_stats = [10 12];
%     s.freqband_perm = [10 12];
    s.timeband_stats = [.6]; s.timeband_perm = [.6];
    s.tf_label_stats = 'Default'; s.tf_label_perm = 'Default';
    [tf_avail] = get_freqband_timeband(s.perm_mode,opts_exclude); s.tf_avail = tf_avail;
    

    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 1;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 1;                % 0=Abs(Diff); 1=Diff;
        opts_pls.perm2pls_split_plusminus = 0;             % 0-Return everything; 1-return pos+negative in separate columns; 2-return only positive (pls.*Ramp(Diff)); 3-return only negative (pls.*Ramp(-Diff)); 4-return only significant cells (doesn't really have a use now)
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 1;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        opts_pls.do_diff_percent = 0;
        opts_pls.do_abs_diff = 0;            % Take absolute value after doing diff.
    opts_pls.target_pls_format = 0; % Convert pls to match this format!
    opts_pls.collapse_pls_to_days = 0;
    opts_pls.spectrogram2spectra_timeslice = 0;   % If working with a spectrogram, take a slice at time given by timeband_stats.
    opts_pls.spectrogram2ts_freqslice = 0;

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
    s.sp_threshold = 40;
    
    % Group options
    s.do_group_collapse_pls2days = 0;
    
    % Groupmode
    s.groupmode = 1;   % 0-Use default grouping (all pairs, enumerate over ctgs);
                 % 1-Group based on permutation test output for diff.
                 % 2-Group based on permutation test output all.
    s.swap_mode = 0;

    s.group_do_merge = 0;
        s.grouppmerge_do_percent = 1;

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
        (opts_pls.permdat2pls == 1 && opts_pls.do_diff == 1 && opts_pls.do_abs_diff == 0);
    
    s.do_custom_colourmap = s.doing_cat_vs_dog;
    s.do_custom_colourmap = false;


% % Plotting options
    paperfig_mode = 1;
    s.opts_PM3Dcs.paperfig_mode=paperfig_mode;
    s.opts_PM3Dcs.stats_mode = 0;
    s.opts_PM3Dcs.do_subplots = 0;
        s.opts_PM3Dcs.max_subplots_per_fig = 16;
    s.opts_PM3Dsp.paperfig_mode=paperfig_mode;
    s.opts_PM3Dsp.symmetric_axes = s.doing_cat_vs_dog;
    s.opts_PM3Dsp.uniform_zaxis_allgroups = 1;           % Makes z-axis the same for all groups plotted
    s.opts_PM3Dsp.do_subplots = 1;
        s.opts_PM3Dsp.max_subplots_per_fig = 16;
    s.opts_PM3Dsp.show_range_stats = 1;
    s.opts_PM3Dsp.show_range_perm = 1;
    s.opts_PSC.paperfig_mode = paperfig_mode;
    s.opts_PSC.remove_dependent = 0;
    s.opts_PSC.do_stats = 1;
    s.opts_PSC.do_diagtest = 0;
    s.opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    s.opts_PSC.opts_diagtest.threshold = 1.0; 
    s.opts_PDays.do_subplots = 1;
        s.opts_PDays.max_subplots_per_fig = 16;
        
    do_animalL = 0;
    do_animalO = 1;
    
    do_positives = 0;
    do_negatives = 1;
    
    
    % Animal L
    if do_animalL
        opts_exclude.excludeL = 0;
        opts_exclude.excludeO = 1; 
        [tf_avail] = get_freqband_timeband(s.perm_mode,opts_exclude); s.tf_avail = tf_avail;

        % Positive bias (towards Cat/Goc)
        if do_positives
            i=1; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                 s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
            
            i=3; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                 s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
            
            i=4; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                 s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
        end

        % Negative bias (towards Dog/Tad)
        if do_negatives
            i=2; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                 s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
            
            i=5; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                 s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
            
            i=6; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                 s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
        end
    end
        
    % Animal O, 
    if do_animalO
        opts_exclude.excludeL = 1;
        opts_exclude.excludeO = 0; 
        [tf_avail] = get_freqband_timeband(s.perm_mode,opts_exclude); s.tf_avail = tf_avail;

        % Positive bias (towards Cat/Goc; modes 1, 2, 8, 11)
        if do_positives
                i=2;  s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                      s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
                
                i=6;  s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                      s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
                
                i=7;  s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                      s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
                
                i=8;  s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                      s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
        end

        % Negative bias (towards Dog/Tad; modes 3, 9, 10)
        if do_negatives
            i=1; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                 s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
            
            i=5; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
                 s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
            [wrkspc_buffer, out] = Fg_5_01a_unit_bias_vs_FFC_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm); title(s.tf_label_perm);
        end
    end

    
%% Figure 5_01c - Analysis of FFC/power biases grouped by single unit biases
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg5_01c';

% Setup params 
    clear group group0

    data_mode = 45;
    switch data_mode
        case 22         % FFC
            s.sfc_mode =  22.401811101;
            s.perm_mode = s.sfc_mode;
%             s.sfc_mode =  22.401511111;     % Partial FFC
%             s.perm_mode = s.sfc_mode;
%             s.sfc_mode =  22.451411103;     % FFC boundary vs non-boundary
              s.perm_mode = 52.700001001;
        case 23         % FFC spectrogram
            s.sfc_mode =  23.4018111011;
            s.perm_mode = s.sfc_mode;
            s.perm_mode = 52.700001001;
        case 41         % PSD
            s.sfc_mode =  41.601811101;
            s.perm_mode = 41.601811101;
            s.perm_mode = 52.700001001;
        case 45         % PSD spectrogram
            s.sfc_mode =  45.6018111011;
            s.perm_mode = s.sfc_mode;
            s.perm_mode = 52.700001001;
        case 52.2         % Units
            s.sfc_mode  =  52.700201001;
            s.perm_mode =  52.700001001;
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
    
%     i=1; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
%     i=1; s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
    

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

    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 3;      % 0-Either; 1-Both; 2-Positive; 3-Negative
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
    
    % Groupmode
    s.groupmode = 1.25;   % 0-Use default grouping (all pairs, enumerate over ctgs);
                 % 1-Group based on permutation test output for diff.
                 % 2-Group based on permutation test output all.
    s.swap_mode = 0;

    s.group_do_merge = 0;
        s.grouppmerge_do_percent = 0;

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
        (opts_pls.permdat2pls == 1 && opts_pls.do_diff == 1 && opts_pls.do_abs_diff == 0);
    
    s.do_custom_colourmap = s.doing_cat_vs_dog;
    s.do_custom_colourmap = false;


% % Plotting options
    paperfig_mode = 1;
    s.opts_PM3Dcs.paperfig_mode=paperfig_mode;
    s.opts_PM3Dcs.stats_mode = 0;
    s.opts_PM3Dcs.do_subplots = 0;
        s.opts_PM3Dcs.max_subplots_per_fig = 16;
    s.opts_PM3Dsp.paperfig_mode=paperfig_mode;
    s.opts_PM3Dsp.symmetric_axes = s.doing_cat_vs_dog;
    s.opts_PM3Dsp.uniform_zaxis_allgroups = 0;           % Makes z-axis the same for all groups plotted
    s.opts_PM3Dsp.do_subplots = 1;
        s.opts_PM3Dsp.max_subplots_per_fig = 16;
    s.opts_PM3Dsp.show_range_stats = 1;
    s.opts_PM3Dsp.show_range_perm = 1;
    s.opts_PSC.paperfig_mode = paperfig_mode;
    s.opts_PSC.remove_dependent = 0;
    s.opts_PSC.do_stats = 1;
    s.opts_PSC.do_diagtest = 0;
    s.opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    s.opts_PSC.opts_diagtest.threshold = 1.0; 
    s.opts_PDays.do_subplots = 1;
        s.opts_PDays.max_subplots_per_fig = 16;


% Basic settings        
second_difference_mode = 4;
prediff_subplots = 1;

% Animal L
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 

opts_perm.split_plusminus = 2;      % 0-Either; 1-Both; 2-Positive; 3-Negative
[wrkspc_buffer, out] = Fg_5_01c_FFC_bias_grouped_by_unit_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm,second_difference_mode,prediff_subplots);

opts_perm.split_plusminus = 3;      % 0-Either; 1-Both; 2-Positive; 3-Negative
[wrkspc_buffer, out] = Fg_5_01c_FFC_bias_grouped_by_unit_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm,second_difference_mode,prediff_subplots);
    

% Animal O
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 

opts_perm.split_plusminus = 2;      % 0-Either; 1-Both; 2-Positive; 3-Negative
[wrkspc_buffer, out] = Fg_5_01c_FFC_bias_grouped_by_unit_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm,second_difference_mode,prediff_subplots);

opts_perm.split_plusminus = 3;      % 0-Either; 1-Both; 2-Positive; 3-Negative
[wrkspc_buffer, out] = Fg_5_01c_FFC_bias_grouped_by_unit_bias(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm,second_difference_mode,prediff_subplots);
    
%% Figure 5_02a - Analyze raw average FFC across cue, fixation, sample, delay
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg5_02a';

% Setup params 
    clear group group0

    data_mode = 23;
    switch data_mode
        case 23         % PSD spectrogram
            sfc_mode1 = 23.4013111013;          % Spectrogram
            sfc_mode2 =  22.401311101;          % PSD
            
            s.sfc_mode =  sfc_mode1;
            s.perm_mode = sfc_mode1;
        case 45         % PSD spectrogram
            sfc_mode1 = 45.6013111013;          % Spectrogram
            sfc_mode2 =  41.601311101;          % PSD
            
            s.sfc_mode =  sfc_mode1;
            s.perm_mode = sfc_mode1;
    end



    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 0;
    opts_exclude.excludeL = 1;
    opts_exclude.excludeO = 0; 
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
    
%     i=1; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
%     i=1; s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
    

    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
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
    
    % Groupmode
    s.groupmode = 0;   % 0-Use default grouping (all pairs, enumerate over ctgs);
                 % 1-Group based on permutation test output for diff.
                 % 2-Group based on permutation test output all.
    s.swap_mode = 0;

    s.group_do_merge = 0;
        s.grouppmerge_do_percent = 1;

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
        (opts_pls.permdat2pls == 1 && opts_pls.do_diff == 1 && opts_pls.do_abs_diff == 0);
    
    s.do_custom_colourmap = s.doing_cat_vs_dog;
    s.do_custom_colourmap = false;


% % Plotting options
    paperfig_mode = 1;
    s.opts_PM3Dcs.paperfig_mode=paperfig_mode;
    s.opts_PM3Dcs.stats_mode = 0;
    s.opts_PM3Dcs.do_subplots = 0;
        s.opts_PM3Dcs.max_subplots_per_fig = 16;
    s.opts_PM3Dsp.paperfig_mode=paperfig_mode;
    s.opts_PM3Dsp.symmetric_axes = s.doing_cat_vs_dog;
    s.opts_PM3Dsp.uniform_zaxis_allgroups = 1;           % Makes z-axis the same for all groups plotted
    s.opts_PM3Dsp.do_subplots = 1;
        s.opts_PM3Dsp.max_subplots_per_fig = 16;
    s.opts_PM3Dsp.show_range_stats = 0;
    s.opts_PM3Dsp.show_range_perm = 1;
    s.opts_PSC.paperfig_mode = paperfig_mode;
    s.opts_PSC.remove_dependent = 0;
    s.opts_PSC.do_stats = 1;
    s.opts_PSC.do_diagtest = 0;
    s.opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    s.opts_PSC.opts_diagtest.threshold = 1.0; 
    s.opts_PDays.do_subplots = 1;
        s.opts_PDays.max_subplots_per_fig = 16;


        
average_all_ctgs = 1;       % 1-Shows average response across all Scheme A and B trials (ctgsetli 9 and 10 averaged together)
                            % 0-Shows ctgs 1, 2, 3, and 4.
        

% Produces spectrogram plot of average across all trials
s.sfc_mode =  sfc_mode1;
s.perm_mode = sfc_mode1;
s.curr_stage_sfc = 4; s.curr_stage_sp = 4;


[wrkspc_buffer, out] = Fg_5_02a_raw_alltrials_CFSD(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm,average_all_ctgs);


% Plot individual stages
i=0;


% Fixation
s.sfc_mode =  sfc_mode1;
s.perm_mode = sfc_mode1;
s.curr_stage_sfc = 4; s.curr_stage_sp = 4;
opts_pls.spectrogram2spectra_timeslice = 1;
s.timeband_stats = [-0.75];

s.plot_on_spect = 0;
[wrkspc_buffer, out] = Fg_5_02a_raw_alltrials_CFSD(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm,average_all_ctgs);
i=i+1; group(i) = out.group; group(i).legend = 'Cue';

% Cue
s.sfc_mode =  sfc_mode1;
s.perm_mode = sfc_mode1;
s.curr_stage_sfc = 4; s.curr_stage_sp = 4;
opts_pls.spectrogram2spectra_timeslice = 1;
s.timeband_stats = [-0.25];

s.plot_on_spect = 0;
[wrkspc_buffer, out] = Fg_5_02a_raw_alltrials_CFSD(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm,average_all_ctgs);
i=i+1; group(i) = out.group; group(i).legend = 'Fixation';


% Sample
% s.sfc_mode =  sfc_mode2;
% s.perm_mode = sfc_mode2;
% s.curr_stage_sfc = 2; s.curr_stage_sp = 2;
% opts_pls.spectrogram2spectra_timeslice = 0;

s.sfc_mode =  sfc_mode1;
s.perm_mode = sfc_mode1;
s.curr_stage_sfc = 4; s.curr_stage_sp = 4;
opts_pls.spectrogram2spectra_timeslice = 1;
s.timeband_stats = [.3];

s.plot_on_spect = 0;
[wrkspc_buffer, out] = Fg_5_02a_raw_alltrials_CFSD(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm,average_all_ctgs);
i=i+1; group(i) = out.group; group(i).legend = 'Sample';

% Delay
% s.sfc_mode =  sfc_mode2;
% s.perm_mode = sfc_mode2;
% s.curr_stage_sfc = 3; s.curr_stage_sp = 3;
% opts_pls.spectrogram2spectra_timeslice = 0;

s.sfc_mode =  sfc_mode1;
s.perm_mode = sfc_mode1;
s.curr_stage_sfc = 4; s.curr_stage_sp = 4;
opts_pls.spectrogram2spectra_timeslice = 1;
s.timeband_stats = [1.3];


s.plot_on_spect = 0;
[wrkspc_buffer, out] = Fg_5_02a_raw_alltrials_CFSD(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm,average_all_ctgs);
i=i+1; group(i) = out.group; group(i).legend = 'Delay';

figure
opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars', double(length(group) < 10)};
[h1] = plot_matrix3D_custstruct([],group,opts_PM3D,s.opts_PM3Dcs);





%% Figure 5_02b - Analyze trends in FFC time series for SIGNIFICANT electrodes
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg5_02b';

% Setup params 
    clear group group0

    data_mode = 23;
    switch data_mode
        case 23         % FFC spectrogram
            s.sfc_mode =  23.4013111011;
            s.perm_mode = s.sfc_mode;
%             s.perm_mode = 52.700001001;
        case 45         % PSD spectrogram
            s.sfc_mode =  45.6013111010;
            s.perm_mode = s.sfc_mode;
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
    opts_exclude.excludeL = 1;
    opts_exclude.excludeO = 0; 
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
    
    i=10; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
    i=10; s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
    

    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
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
    opts_pls.spectrogram2ts_freqslice = 1;

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
    
    % Groupmode
    s.groupmode = 3;   % 0-Use default grouping (all pairs, enumerate over ctgs);
                 % 1-Group based on permutation test output for diff.
                 % 2-Group based on permutation test output all.
    s.swap_mode = 0;

    s.group_do_merge = 0;
        s.grouppmerge_do_percent = 1;

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
        (opts_pls.permdat2pls == 1 && opts_pls.do_diff == 1 && opts_pls.do_abs_diff == 0);
    
    s.do_custom_colourmap = s.doing_cat_vs_dog;
    s.do_custom_colourmap = false;


% % Plotting options
    paperfig_mode = 1;
    s.opts_PM3Dcs.paperfig_mode=paperfig_mode;
    s.opts_PM3Dcs.stats_mode = 0;
    s.opts_PM3Dcs.do_subplots = 0;
        s.opts_PM3Dcs.max_subplots_per_fig = 16;
    s.opts_PM3Dsp.paperfig_mode=paperfig_mode;
    s.opts_PM3Dsp.symmetric_axes = s.doing_cat_vs_dog;
    s.opts_PM3Dsp.uniform_zaxis_allgroups = 1;           % Makes z-axis the same for all groups plotted
    s.opts_PM3Dsp.do_subplots = 1;
        s.opts_PM3Dsp.max_subplots_per_fig = 16;
    s.opts_PM3Dsp.show_range_stats = 1;
    s.opts_PM3Dsp.show_range_perm = 1;
    s.opts_PSC.paperfig_mode = paperfig_mode;
    s.opts_PSC.remove_dependent = 0;
    s.opts_PSC.do_stats = 1;
    s.opts_PSC.do_diagtest = 0;
    s.opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    s.opts_PSC.opts_diagtest.threshold = 1.0; 
    s.opts_PDays.do_subplots = 1;
        s.opts_PDays.max_subplots_per_fig = 16;
    
        
        
    % Spectrograms
    opts_pls.spectrogram2ts_freqslice = 0;
    [wrkspc_buffer, out] = Fg_5_02b_raw_sigElects(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
        
    % Freqband slices
    opts_pls.spectrogram2ts_freqslice = 1;
    [wrkspc_buffer, out] = Fg_5_02b_raw_sigElects(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);

    


    % Plot timeband slices
    
    
    
    chosen_sch_sig = 3;            % For plotting responses of significant
    chosen_sch_nonsig = 6;         % For plotting responses of non-sig cells
    
    overlay_non_significant_electrode_responses = 1;
    
    opts_pls.spectrogram2spectra_timeslice = 1;   % If working with a spectrogram, take a slice at time given by timeband_stats.
    opts_pls.spectrogram2ts_freqslice = 0;
    s.plot_on_spect = 0;
    i=0; clear group

    s.timeband_stats = [-.75];
    [wrkspc_buffer, out] = Fg_5_02b_raw_sigElects(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
    i=i+1; group(i) = out.group(chosen_sch_sig); group(i).legend = ['Cue ' group(i).legend];
           group2(i) = out.group(chosen_sch_nonsig); group2(i).legend = ['Cue ' group(i).legend];
    
    s.timeband_stats = [-.25];
    [wrkspc_buffer, out] = Fg_5_02b_raw_sigElects(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
    i=i+1; group(i) = out.group(chosen_sch_sig); group(i).legend = ['Fixation ' group(i).legend];
           group2(i) = out.group(chosen_sch_nonsig); group2(i).legend = ['Fixation ' group(i).legend];

    s.timeband_stats = [.3];
    [wrkspc_buffer, out] = Fg_5_02b_raw_sigElects(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
    i=i+1; group(i) = out.group(chosen_sch_sig); group(i).legend = ['Sample ' group(i).legend];
           group2(i) = out.group(chosen_sch_nonsig); group2(i).legend = ['Sample ' group(i).legend];

    s.timeband_stats = [1.3];
    [wrkspc_buffer, out] = Fg_5_02b_raw_sigElects(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm);
    i=i+1; group(i) = out.group(chosen_sch_sig); group(i).legend = ['Delay ' group(i).legend];
           group2(i) = out.group(chosen_sch_nonsig); group2(i).legend = ['Delay ' group(i).legend];


    figure
    opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars', 1};
    [h1] = plot_matrix3D_custstruct([],group,opts_PM3D,s.opts_PM3Dcs);
    
    
    if overlay_non_significant_electrode_responses
        hold on;
        opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars', 0};
        [h1] = plot_matrix3D_custstruct([],group2,opts_PM3D,s.opts_PM3Dcs);
    end


