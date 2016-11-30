
% New plotting function that uses generalized code
%   See plots_generalized for template or funcs_generalized for useful
%   functions

%%
error('This is meant to run in cell mode. Dont F5!');


%% All figures

% Clear
clearvars -except wrkspc_buffer fv fv3

addpath('figs_batch_4_0_supporting');
addpath(genpath('funcs_plots_preferred'));
addpath(genpath('funcs_generalized'));

% Load defaults
fv3 = figs_batch_defaults_3; 

% Create wrkspc_buffer structure if not already present
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end



%% Figure 4_00 - Alpha and beta peaks raw


clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg4_00';

% Set up pls fv (popts)
% Unit exclusion
opts_exclude = Opts_Exclude;
opts_exclude.exclude_clipping = 1;
opts_exclude.exclude_60 = 0;
opts_exclude.exclude_nans = 1;
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 0; 
opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC

% Stage selection
curr_stage_sfc = 2;
curr_stage_sp = 2;
[wrkspc_buffer, group_s] = Fg_4_00_raw_pls(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,opts_exclude,'Sample');

% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 3;
[wrkspc_buffer, group_d] = Fg_4_00_raw_pls(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,opts_exclude,'Delay');

inds = 9;
figure;
[h1] = plot_matrix3D_custstruct([],[group_s(inds),group_d(inds)],opts_PM3D,opts_PM3Dcs);


%% Figure 4_01 - Beta sens split vs mean alpha pls

clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg4_01';

plot_on_spect = 1;
plot_on_bargraph = 1;
plot_on_func = 0;

do_units = 0;

% SFC mode
if do_units==0
    sfc_mode = 22.401311101;
    perm_mode = 22.401311101;
% sfc_mode = 2.3015111;
% perm_mode = 2.3015111;
% perm_mode = 5;
% sfc_mode = 41.6514111;
% perm_mode = 41.6514111;
elseif do_units == 1
    sfc_mode = 52.7502010;
    perm_mode = 52.7500010;
elseif do_units == 2
    sfc_mode = 52.7500010;
    perm_mode = 52.7500010;
end

% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 3;

if do_units == 1
    curr_stage_sfc = 4;
end

    % Pls switches
    freqband_stats = [16 20];
    freqband_stats_perm = [16 20];
%     freqband_stats = [1 10];
%     freqband_stats_perm = [1 10];
%     freqband_stats = [10 12];
%     freqband_stats_perm = [10 12];
%     freqband_stats = [20 30];
%     freqband_stats_perm = [20 30];

% Set up pls fv (popts)
% Unit exclusion
opts_exclude = Opts_Exclude;
opts_exclude.exclude_clipping = 1;
opts_exclude.exclude_60 = 0;
opts_exclude.exclude_nans = 1;
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 
opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC

% Permutation test options
opts_perm = Opts_Perm;
opts_perm.do_bh0 = 1;
opts_perm.do_phi = 0;
opts_perm.split_plusminus = 1;
opts_perm.alpha0 = 0.05;                 % *********** Setting this huge to get all traces!! *********
opts_perm.alpha_bh0 = 0.8;              % *********** Setting this huge to get all traces!! *********



% Animal L
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 

groupmode = 3;                  % 1/3 Sch A (pooled / notpooled)
[wrkspc_buffer, group] = Fg_4_01_alpha_pls(wrkspc_buffer,sfc_mode, perm_mode,opts_exclude, opts_perm,freqband_stats,freqband_stats_perm,groupmode)
% groupmode = 4;                  % 1/3 Sch A (pooled / notpooled)
% [wrkspc_buffer, group] = Fg_4_01_alpha_pls(wrkspc_buffer,sfc_mode, perm_mode,opts_exclude, opts_perm,freqband_stats,freqband_stats_perm,groupmode)


% Animal O
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 
groupmode = 5;                  % 2/4 Sch B (pooled / notpooled)
[wrkspc_buffer, group] = Fg_4_01_alpha_pls(wrkspc_buffer,sfc_mode, perm_mode,opts_exclude, opts_perm,freqband_stats,freqband_stats_perm,groupmode)
% groupmode = 6;                  % 2/4 Sch B (pooled / notpooled)
% [wrkspc_buffer, group] = Fg_4_01_alpha_pls(wrkspc_buffer,sfc_mode, perm_mode,opts_exclude, opts_perm,freqband_stats,freqband_stats_perm,groupmode)


%% Figure 4_01b - Beta bias, phi analysis

clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg4_01b';

plot_on_spect = 1;
plot_on_bargraph = 1;
plot_on_func = 0;

do_units = 0;

% SFC mode
if do_units==0
    sfc_mode = 22.401311101;
    perm_mode = 22.401311101;
% sfc_mode = 2.3015111;
% perm_mode = 2.3015111;
% perm_mode = 5;
% sfc_mode = 41.6514111;
% perm_mode = 41.6514111;
elseif do_units == 1
    sfc_mode = 52.7502010;
    perm_mode = 52.7500010;
elseif do_units == 2
    sfc_mode = 52.7500010;
    perm_mode = 52.7500010;
end

% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 3;

if do_units == 1
    curr_stage_sfc = 4;
end

    % Pls switches
    freqband_stats = [16 20];
    freqband_stats_perm = [16 20];
%     freqband_stats = [1 10];
%     freqband_stats_perm = [1 10];
    freqband_stats = [9 9];
%     freqband_stats_perm = [10 12];
%     freqband_stats = [20 30];
%     freqband_stats_perm = [20 30];

% Set up pls fv (popts)
% Unit exclusion
opts_exclude = Opts_Exclude;
opts_exclude.exclude_clipping = 1;
opts_exclude.exclude_60 = 0;
opts_exclude.exclude_nans = 1;
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 
opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC

% Permutation test options
opts_perm = Opts_Perm;
opts_perm.do_bh0 = 1;
opts_perm.do_phi = 0;
opts_perm.split_plusminus = 1;
opts_perm.alpha0 = 0.05;                 % *********** Setting this huge to get all traces!! *********
opts_perm.alpha_bh0 = 0.5;              % *********** Setting this huge to get all traces!! *********



% Animal L
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 
groupmode = 1;                  % 1/3 Sch A (pooled / notpooled)
[wrkspc_buffer, group] = Fg_4_01b_beta_bias_phi(wrkspc_buffer,sfc_mode, perm_mode,opts_exclude, opts_perm,freqband_stats,freqband_stats_perm,groupmode);


% Animal O
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 
groupmode = 2;                  % 2/4 Sch B (pooled / notpooled)
[wrkspc_buffer, group] = Fg_4_01b_beta_bias_phi(wrkspc_buffer,sfc_mode, perm_mode,opts_exclude, opts_perm,freqband_stats,freqband_stats_perm,groupmode);


%% Figure 4_02 - Beta sens split vs alpha z-score
% This, unfortunately, doesn't give any good results. But its the wrong
% idea. Instead, assume all "fat" electrodes are doing alpha and "fat beta"
% electrodes are simply a byproduct of this. 
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg4_02';

plot_on_spect = 1;
plot_on_bargraph = 1;
plot_on_func = 0;

do_units = 0;
perm_mode = 22.401811101;

% SFC mode
if do_units==0
    sfc_mode = 22.401811101;    
% sfc_mode = 2.3015111;
% perm_mode = 2.3015111;
% perm_mode = 5;
% sfc_mode = 41.6514111;
% perm_mode = 41.6514111;
elseif do_units == 1
    sfc_mode = 52.7002010;
elseif do_units == 2
    sfc_mode = 52.7000010;
end

% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 3;

if do_units == 1
    curr_stage_sfc = 4;
end

    % Pls switches
    freqband_stats = [16 20];
    freqband_stats_perm = [16 20];
%     freqband_stats = [1 10];
%     freqband_stats_perm = [1 10];
    freqband_stats = [10 12];
%     freqband_stats_perm = [10 12];
%     freqband_stats = [20 30];
%     freqband_stats_perm = [20 30];

% Set up pls fv (popts)
% Unit exclusion
opts_exclude = Opts_Exclude;
opts_exclude.exclude_clipping = 1;
opts_exclude.exclude_60 = 0;
opts_exclude.exclude_nans = 1;
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 
opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC

% Permutation test options
opts_perm = Opts_Perm;
opts_perm.do_bh0 = 0;
opts_perm.do_phi = 0;
opts_perm.split_plusminus = 1;
opts_perm.alpha0 = 0.05;
opts_perm.alpha_bh0 = 0.2;


%
% Animal L
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 
groupmode = 3;                  % 0-SchA&B; 1-SchA; 2-SchB
[wrkspc_buffer, group] = Fg_4_02_alpha_perm(wrkspc_buffer,sfc_mode, perm_mode,opts_exclude, opts_perm,freqband_stats,freqband_stats_perm,groupmode,curr_stage_sfc,curr_stage_sp)


% Animal O
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 
groupmode = 4;                  % 0-SchA&B; 1-SchA; 2-SchB
[wrkspc_buffer, group] = Fg_4_02_alpha_perm(wrkspc_buffer,sfc_mode, perm_mode,opts_exclude, opts_perm,freqband_stats,freqband_stats_perm,groupmode,curr_stage_sfc,curr_stage_sp)



%% Figure 4_03 - Beta vs alpha scatterplot
% This, unfortunately, doesn't give any interpretable results
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg4_03';

plot_on_spect = 1;
plot_on_bargraph = 1;
plot_on_func = 0;

do_units = 0;

% SFC mode
if do_units==0
    sfc_mode = 22.401811101;
    perm_mode = 22.401811101;
% sfc_mode = 2.3015111;
% perm_mode = 2.3015111;
% perm_mode = 5;
% sfc_mode = 41.6514111;
% perm_mode = 41.6514111;
elseif do_units == 1
    sfc_mode = 52.7502010;
    perm_mode = 52.7500010;
elseif do_units == 2
    sfc_mode = 52.7500010;
    perm_mode = 52.7500010;
end

% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 3;

if do_units == 1
    curr_stage_sfc = 4;
end

    % Pls switches
    freqband_stats2 = [16 20];
    freqband_stats_perm2 = [16 20];
%     freqband_stats = [1 10];
%     freqband_stats_perm = [1 10];
    freqband_stats1 = [10 12];
    freqband_stats_perm1 = [10 12];
%     freqband_stats = [20 30];
%     freqband_stats_perm = [20 30];

% Set up pls fv (popts)
% Unit exclusion
opts_exclude = Opts_Exclude;
opts_exclude.exclude_clipping = 1;
opts_exclude.exclude_60 = 0;
opts_exclude.exclude_nans = 1;
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 
opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC

% Permutation test options
opts_perm = Opts_Perm;
opts_perm.do_bh0 = 1;
opts_perm.do_phi = 0;
opts_perm.split_plusminus = 1;
opts_perm.alpha0 = 0.05;
opts_perm.alpha_bh0 = 0.2;


% Animal L - Cat / Dog
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 
groupmode = 3;                  % 0-SchA&B; 1-SchA; 2-SchB
[wrkspc_buffer, group] = Fg_4_03_scatterplot_beta_alpha(wrkspc_buffer,sfc_mode, perm_mode,opts_exclude, opts_perm,freqband_stats1,freqband_stats_perm1,freqband_stats2,freqband_stats_perm2,groupmode)

%
% Animal O - Fat / Thin
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 
groupmode = 4;                  % 0-SchA&B; 1-SchA; 2-SchB
[wrkspc_buffer, group] = Fg_4_03_scatterplot_beta_alpha(wrkspc_buffer,sfc_mode, perm_mode,opts_exclude, opts_perm,freqband_stats1,freqband_stats_perm1,freqband_stats2,freqband_stats_perm2,groupmode)


%% Figure 4_03b - Beta, alpha, unit histograms
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg4_03b';

do_units = 0;

% SFC mode
if do_units==0
    sfc_mode = 22.401811101;
    perm_mode = 22.401811101;
%     sfc_mode = 41.6014111;
%     perm_mode = 41.6014111;
elseif do_units == 1
    sfc_mode = 52.7002010;
    perm_mode = 52.7000010;
elseif do_units == 2
    sfc_mode = 52.7000010;
    perm_mode = 52.7000010;
elseif do_units == 4        % Spectrogram
    sfc_mode = 23.4014111;
    perm_mode = 23.4014111;
end


% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 3;

if any(do_units == [1,4])
    curr_stage_sfc = 4;
end

    % Pls switches
    freqband_stats = [16 20];
    freqband_stats_perm = [16 20];
%     freqband_stats = [10 12];
%     freqband_stats_perm = [10 12];

    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 0;
    opts_exclude.excludeL = 1;
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC

% Histogram settings
hs.nbins_or_edges = [];
% hs.nbins_or_edges = 20;    % Use 100 for FFC; comment out for units
hs.show_vertical_bar = true;
hs.symmetric_about_x_axis = true;
xldesired = [];
% xldesired  = [-.3 .3];    % For FFC
% xldesired  = [-20 20];    % For units (15 for sample; 20 for delay)

% Animal L - Cat vs Dog
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 
    % Beta
    freqband_stats = [16 20]; freqband_stats_perm = [16 20];
    groupmode=3;
    [wrkspc_buffer, group] = Fg_4_03b_hist(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,opts_exclude,groupmode)
    figure; plot_histplot(group,'plot_type','hist_autogrouped','bar_orientation','vert','hist_settings',hs);
    colormap([get_clist([1:2],'matrix'); 0.2 0.2 0.2]);
    if ~isempty(xldesired); xlim(xldesired); end
    
    freqband_stats = [16 20]; freqband_stats_perm = [16 20];
    groupmode=4;
    [wrkspc_buffer, group] = Fg_4_03b_hist(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,opts_exclude,groupmode)
    figure; plot_histplot(group,'plot_type','hist_autogrouped','bar_orientation','vert','hist_settings',hs);
    colormap([get_clist([3:4],'matrix'); 0.2 0.2 0.2]);
    if ~isempty(xldesired); xlim(xldesired); end
    
%     % Alpha
%     freqband_stats = [10 12]; freqband_stats_perm = [10 12];
%     groupmode=3;
%     [wrkspc_buffer, group] = Fg_4_03b_hist(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,opts_exclude,groupmode)
%     figure; plot_histplot(group,'plot_type','hist_autogrouped','bar_orientation','vert','hist_settings',hs);
%     colormap([get_clist([1:2],'matrix'); 0.2 0.2 0.2]);
%     if ~isempty(xldesired); xlim(xldesired); end
%     
%     freqband_stats = [10 12]; freqband_stats_perm = [10 12];
%     groupmode=4;
%     [wrkspc_buffer, group] = Fg_4_03b_hist(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,opts_exclude,groupmode)
%     figure; plot_histplot(group,'plot_type','hist_autogrouped','bar_orientation','vert','hist_settings',hs);
%     colormap([get_clist([3:4],'matrix'); 0.2 0.2 0.2]);
%     if ~isempty(xldesired); xlim(xldesired); end
    
% Animal O - Fat vs Thin
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 
    % Beta
    freqband_stats = [16 20]; freqband_stats_perm = [16 20];
    groupmode=3;
    [wrkspc_buffer, group] = Fg_4_03b_hist(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,opts_exclude,groupmode)
    figure; plot_histplot(group,'plot_type','hist_autogrouped','bar_orientation','vert','hist_settings',hs);
    colormap([get_clist([1:2],'matrix'); 0.2 0.2 0.2]);
    if ~isempty(xldesired); xlim(xldesired); end
    
    freqband_stats = [16 20]; freqband_stats_perm = [16 20];
    groupmode=4;
    [wrkspc_buffer, group] = Fg_4_03b_hist(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,opts_exclude,groupmode)
    figure; plot_histplot(group,'plot_type','hist_autogrouped','bar_orientation','vert','hist_settings',hs);
    colormap([get_clist([3:4],'matrix'); 0.2 0.2 0.2]);
    if ~isempty(xldesired); xlim(xldesired); end
    
%     % Alpha
%     freqband_stats = [10 12]; freqband_stats_perm = [10 12];
%     groupmode=3;
%     [wrkspc_buffer, group] = Fg_4_03b_hist(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,opts_exclude,groupmode)
%     figure; plot_histplot(group,'plot_type','hist_autogrouped','bar_orientation','vert','hist_settings',hs);
%     colormap([get_clist([1:2],'matrix'); 0.2 0.2 0.2]);
%     if ~isempty(xldesired); xlim(xldesired); end
%     
%     freqband_stats = [10 12]; freqband_stats_perm = [10 12];
%     groupmode=4;
%     [wrkspc_buffer, group] = Fg_4_03b_hist(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,opts_exclude,groupmode)
%     figure; plot_histplot(group,'plot_type','hist_autogrouped','bar_orientation','vert','hist_settings',hs);
%     colormap([get_clist([3:4],'matrix'); 0.2 0.2 0.2]);
%     if ~isempty(xldesired); xlim(xldesired); end


%% Figure 4_04 - FFC vs Units for Fat/Thin


clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg4_04';

% SFC mode
sfc_mode = 22.4014111;
% perm_mode = 22.4014111;
% sfc_mode =  2.3515111;
% perm_mode = 2.3015111;
% perm_mode = 5;
% sfc_mode = 41.6514111;
% perm_mode = 41.6514111;
% sfc_mode = 52.7002010;
perm_mode = 52.7000010;


% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 3;

    % Pls switches
    freqband_stats = [16 20];
    freqband_stats_perm = [16 20];
%     freqband_stats = [1 10];
%     freqband_stats_perm = [1 10];
%     freqband_stats = [10 12];
%     freqband_stats_perm = [9 11];
%     freqband_stats_perm = [8 10];
%     freqband_stats_perm = [10 12];
%     freqband_stats = [30 40];
%     freqband_stats_perm = [30 40];
%     freqband_stats_perm = [80 100];

    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 1;
    opts_exclude.excludeL = 1;
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC

    warning('To do - implement Buschmann figure.');
    warning('Add horizontal bar in perm figures to indicate value expected by chance');
    warning('Remove alpha vertical line');
    warning('Redo spectrograms without showing T stage');

    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.perm2pls_split_plusminus = 0;             % 0-Return everything; 1-return pos+negative in separate columns; 2-return only positive (pls.*Ramp(Diff)); 3-return only negative (pls.*Ramp(-Diff)); 4-return only significant cells (doesn't really have a use now)
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))

    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 0;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 1;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dcs.stats_mode = 0;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    
    warning('To do - compare sch sensitive electrodes ability to attenuate irrelevant ctg response to other electrodes.');

% Plot switches
plot_on_spect = 1;
plot_on_spectrogram = 0;
plot_on_scatter = 0;
plot_on_bargraph = 0;


group_do_merge = 0;
    grouppmerge_do_percent = 1;


        
% Load pls
[wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
pls = out_pls.pls;
pls_stats = out_pls.pls_stats;
abscissa = out_pls.abscissa;
abscissa2 = out_pls.abscissa2;
bad_any = out_pls.bad_any;
funames1D = out_pls.funames1D;
mypairs = out_pls.mypairs;
group0 = out_pls.group;
clear out_pls

% pls = pls * -1;
% pls_stats = pls_stats * -1;


% Build group
groupmode = 1;

% Use default grouping (all pairs, enumerate over ctgs)
switch groupmode
    case 0
        group = group0;
        sp = ~bad_any(:);
    case 1
    % Load sp's
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm);
    sp = out_perm.sig_cells;


    % Load bads perm
    [bad_any_perm] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);

    % Map sp's as needed
    [sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any);

    % Create group template
    mycrit = [2*ones(1,size(sp,2))];
    grt = group0(1);
    grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;

    % Run a simple test
    clear group
    i=0;
    %i=i+1; group(i)=grt; group(i).criteria(1:2)=[1 0]; group(i).ctgs=1;   % Ctg1-2 deciders
    %i=i+1; group(i)=grt; group(i).criteria(1:2)=[0 0]; group(i).ctgs=1;   % Ctg1-2 non-deciders
    
%     i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 0 1 0]; group(i).ctgs=2;   % Ctg1-2 deciders
%     i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 0 0 0]; group(i).ctgs=2;   % Ctg1-2 non-deciders
    i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 0 0 1]; group(i).ctgs=2;   % Ctg1-2 deciders
    i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 0 0 0]; group(i).ctgs=2;   % Ctg1-2 non-deciders

end

% Load data into groups
group = group.query_legend(group0);
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);

if group_do_merge
    % Flip if ctgsetli mode = 4 or 5.
    temp=1:length(group);
    N=length(group);
    temp=flipud(reshape(1:N,2,floor(N/2)));
    temp=temp(:);
    group = grouppairs_merge(group(temp),opts_pls.perm2pls_dophi,grouppmerge_do_percent);
end

% Test Plot groups, at last
opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
if plot_on_spect
    
    
    inds = 1:length(group);
%     inds = [1:4];
%     inds = [5:8];
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    
end


% Test spectrogram
if plot_on_spectrogram
    plot_spectrogram_custstruct(group(:),opts_PM3Dsp);
    
    
end


% Test bargraph

if plot_on_bargraph
    fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
end


% Test Plot scattergroups

if plot_on_scatter
    figure; % subplot(221); 
    ind1 = 1;
    ind2 = 3;
    x = pls_stats(~bad_any,group(ind1).ctgs);
    y = pls_stats(~bad_any,group(ind2).ctgs);
    plott_fit(x,y,'k.');
    hold on; plot_scattergroups(x,y,[group(ind1),group(ind2)],~bad_any);
    xlabel(['FFC ' group(ind1).legend]); ylabel(['FFC ' group(ind2).legend]); % Might need fixing
    
    
    fprintf(['L1 Slope SchA Rel Vs Irrel= ' num2str(std(pls_stats(~bad_any,3))/std(pls_stats(~bad_any,1))) '\n']);
    fprintf(['L1 Slope SchB Rel Vs Irrel= ' num2str(std(pls_stats(~bad_any,4))/std(pls_stats(~bad_any,2))) '\n']);
    
    
end

do_pca = 0;
if do_pca
   [coeff,score,latent] = princomp([x(:),y(:)]);
    ax1 = coeff(:,1);% *sqrt(latent(1));
    ax2 = coeff(:,2);% *sqrt(latent(2));
    hold on; plot([0 coeff(1,1)],[0 coeff(2,1)],'r','LineWidth',2);
    hold on; plot([0 coeff(1,2)],[0 coeff(2,2)],'r','LineWidth',2); 
end

%% Figure 4_05 - FFC perm split_plusminus
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg4_05';

do_units = 0;

% SFC mode
if do_units==0
    sfc_mode = 22.401811101;
    perm_mode = 22.401811101;
%     sfc_mode = 41.601311101;
%     perm_mode = 41.601311101;
elseif do_units == 1
    sfc_mode = 52.7002010;
    perm_mode = 52.7000010;
elseif do_units == 2
    sfc_mode = 52.7000010;
    perm_mode = 52.7000010;
elseif do_units == 4        % Spectrogram
    sfc_mode = 23.4014111;
    perm_mode = 23.4014111;
end


% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 3;

if any(do_units == [1,4])
    curr_stage_sfc = 4;
end

plot_on_spect = 1;
plot_on_spectrogram = 0;
plot_on_scatter = 0;
plot_on_bargraph = 1;

switch do_units
    case 2; plot_on_spect = 0;
    case 4;
        plot_on_spectrogram = 1;
        plot_on_spect = 0;
        plot_on_bargraph = 0;
end


    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 1;
    opts_exclude.excludeL = 1;
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % % % % % Monkey L % % % % % 
    opts_exclude.excludeL = 0;
    opts_exclude.excludeO = 1; 
    
    % No splitting
    groupmode=0;
    [wrkspc_buffer, group] = Fg_4_05_FFC_perm_split_plusminus(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,groupmode);
    
    % With splitting
    groupmode=1;
    [wrkspc_buffer, group_s] = Fg_4_05_FFC_perm_split_plusminus(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,groupmode);
    
    if plot_on_spect
        yl = [];
%         yl = [0 120];
        inds = [1]; inds_s = [1,2]; group(inds).ylims_desired = [];
        figure; [h1] = plot_matrix3D_custstruct([],[group(inds) group_s(inds_s)],opts_PM3D,opts_PM3Dcs);
        xl = xlim; hold on; plot([xl(1) xl(2)], [1 1].*group(inds).numcells*0.05/2,'k:','LineWidth',2);   % Add horizontal line for expect
        
        inds = [2]; inds_s = [3,4]; group(inds).ylims_desired = yl;
        figure; [h1] = plot_matrix3D_custstruct([],[group(inds) group_s(inds_s)],opts_PM3D,opts_PM3Dcs);
        xl = xlim; hold on; plot([xl(1) xl(2)], [1 1].*group(inds).numcells*0.05/2,'k:','LineWidth',2);   % Add horizontal line for expect
    end

    
    % Test spectrogram
    if plot_on_spectrogram
        plot_spectrogram_custstruct(group(plot_inds),opts_PM3Dsp);
    end

    
    % % % % % Monkey O % % % % % 
    opts_exclude.excludeL = 1;
    opts_exclude.excludeO = 0; 
    
    % No splitting
    groupmode=0;
    [wrkspc_buffer, group] = Fg_4_05_FFC_perm_split_plusminus(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,groupmode);
    
    % With splitting
    groupmode=1;
    [wrkspc_buffer, group_s] = Fg_4_05_FFC_perm_split_plusminus(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,groupmode);
    
    if plot_on_spect
        yl = [];
%         yl = [0 400];
        inds = [1]; inds_s = [1,2]; group(inds).ylims_desired = [];
        figure; [h1] = plot_matrix3D_custstruct([],[group(inds) group_s(inds_s)],opts_PM3D,opts_PM3Dcs);
        xl = xlim; hold on; plot([xl(1) xl(2)], [1 1].*group(inds).numcells*0.05/2,'k:','LineWidth',2);   % Add horizontal line for expect
        
        inds = [2]; inds_s = [3,4]; group(inds).ylims_desired = [yl];
        figure; [h1] = plot_matrix3D_custstruct([],[group(inds) group_s(inds_s)],opts_PM3D,opts_PM3Dcs);
        xl = xlim; hold on; plot([xl(1) xl(2)], [1 1].*group(inds).numcells*0.05/2,'k:','LineWidth',2);   % Add horizontal line for expect
    end


    % Test spectrogram
    if plot_on_spectrogram
        plot_spectrogram_custstruct(group(plot_inds),opts_PM3Dsp);
    end
    
%% Figure 4_06 - Beta FFC versus beta PSD correlation

clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg4_06';

% More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.perm2pls_split_plusminus = 0;             % 0-Return everything; 1-return pos+negative in separate columns; 2-return only positive (pls.*Ramp(Diff)); 3-return only negative (pls.*Ramp(-Diff)); 4-return only significant cells (doesn't really have a use now)
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))

                
  
    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 1;
    opts_exclude.excludeL = 0;
    opts_exclude.excludeO = 1; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    plotmode1 = 1;
    plotmode2 = 1;
 
warning('Next round'); pause
% % % % % % % % % % % % % % % % % % % Round 1 - Beta FFC vs beta PSD - Animal L only % % % % % 
% % % % % % % % % % % % % % % % % % Sample stage
% Animal selection
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 

sfc_mode1 = 22.4013111;
sfc_mode2 = 41.6013111;
curr_stage1 = 2;
curr_stage2 = 2;
freqband_stats1 = [16 20];
freqband_stats2 = [16 20];
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);

% % % % % % % % % % % % % % % % % % Delay stage
% Animal selection
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 

sfc_mode1 = 22.4013111;
sfc_mode2 = 41.6013111;
curr_stage1 = 3;
curr_stage2 = 3;
freqband_stats1 = [16 20];
freqband_stats2 = [16 20];
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);



warning('Next round'); pause
% % % % % % % % % % % % % % % % % % % Round 2 - Alpha FFC vs Alpha PSD - Delay stage only % % % % % 
% % % % % % % % % % % % % % % % % % Animal L
% Animal selection
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 

sfc_mode1 = 22.4013111;
sfc_mode2 = 41.6013111;
curr_stage1 = 3;
curr_stage2 = 3;
freqband_stats1 = [10 12];
freqband_stats2 = [10 12];
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);

% % % % % % % % % % % % % % % % % % Animal O
% Animal selection
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 

sfc_mode1 = 22.4013111;
sfc_mode2 = 41.6013111;
curr_stage1 = 3;
curr_stage2 = 3;
freqband_stats1 = [10 12];
freqband_stats2 = [10 12];
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);

warning('Next round'); pause
% % % % % % % % % % % % % % % % % % % Round 3 - Alpha vs beta PSD - Animal L delay only % % % % % 
% % % % % % % % % % % % % % % % % % Delay stage, animal L
% Animal selection
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 

sfc_mode1 = 41.6013111;
sfc_mode2 = 41.6013111;
curr_stage1 = 3;
curr_stage2 = 3;
freqband_stats1 = [10 12];
freqband_stats2 = [16 20];
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);

warning('Next round'); pause
% % % % % % % % % % % % % % % % % % % Round 4 - Alpha vs beta FFC % % % % % 
% % % % % % % % % % % % % % % % % % Animal L
% Animal selection
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 

sfc_mode1 = 22.4013111;
sfc_mode2 = 22.4013111;
curr_stage1 = 3;
curr_stage2 = 3;
freqband_stats1 = [10 12];
freqband_stats2 = [16 20];
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);

% % % % % % % % % % % % % % % % % % Animal O
% Animal selection
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 

sfc_mode1 = 22.4013111;
sfc_mode2 = 22.4013111;
curr_stage1 = 3;
curr_stage2 = 3;
freqband_stats1 = [10 12];
freqband_stats2 = [16 20];
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);


warning('Next round'); pause
% % % % % % % % % % % % % % % % % % % Round 5 - Alpha PSD vs Beta FFC % % % % % 
% % % % % % % % % % % % % % % % % % Animal O, Sample Stage, Alpha PSD vs Beta FFC
% Animal selection
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0;

sfc_mode1 = 22.4013111;
sfc_mode2 = 41.6013111;
curr_stage1 = 2;
curr_stage2 = 2;
freqband_stats1 = [16 20];
freqband_stats2 = [10 12];
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);

% % % % % % % % % % % % % % Animal L, Delay Stage, Alpha PSD vs Beta FFC
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1;

sfc_mode1 = 22.4013111;
sfc_mode2 = 41.6013111;
curr_stage1 = 3;
curr_stage2 = 3;
freqband_stats1 = [16 20];
freqband_stats2 = [10 12];
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);

% % % % % % % % % % % % % % Animal O, Delay Stage, Alpha PSD vs Beta FFC
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0;

sfc_mode1 = 22.4013111;
sfc_mode2 = 41.6013111;
curr_stage1 = 3;
curr_stage2 = 3;
freqband_stats1 = [16 20];
freqband_stats2 = [10 12];
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);


warning('Next round'); pause
% % % % % % % % % % % % % % % % % % Animal O
% Animal selection
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 

sfc_mode1 = 22.4013111;
sfc_mode2 = 22.4013111;
curr_stage1 = 3;
curr_stage2 = 3;
freqband_stats1 = [10 12];
freqband_stats2 = [16 20];
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);

%% Figure 4_07 - FFC Amp vs Phase

clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg4_07';

% More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.perm2pls_split_plusminus = 0;             % 0-Return everything; 1-return pos+negative in separate columns; 2-return only positive (pls.*Ramp(Diff)); 3-return only negative (pls.*Ramp(-Diff)); 4-return only significant cells (doesn't really have a use now)
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))

                
  
    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 1;
    opts_exclude.excludeL = 0;
    opts_exclude.excludeO = 1; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
 
warning('Next round'); pause
% % % % % % % % % % % % % % % % % % % Round 1 - Beta FFC Amp versus Phase % % % % % 
% % % % % % % % % % % % % % % % % % Sample Animal L
% Animal selection
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 

sfc_mode1 = 22.4013111;
sfc_mode2 = 22.4013111;
curr_stage1 = 2;
curr_stage2 = 2;
freqband_stats1 = [16 20];
freqband_stats2 = [16 20];
plotmode1 = 1;
plotmode2 = 4;
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);

warning('Next round'); pause
% % % % % % % % % % % % % % % % % % Delay Animal L
% Animal selection
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 

sfc_mode1 = 22.4018111;
sfc_mode2 = 22.4018111;
curr_stage1 = 3;
curr_stage2 = 3;
freqband_stats1 = [16 20];
freqband_stats2 = [16 20];
plotmode1 = 1;
plotmode2 = 4;
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);

%
% % % % % % % % % % % % % % % % % % Delay Animal O
% Animal selection
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 

sfc_mode1 = 22.4018111;
sfc_mode2 = 22.4018111;
curr_stage1 = 3;
curr_stage2 = 3;
freqband_stats1 = [16 20];
freqband_stats2 = [16 20];
plotmode1 = 1;
plotmode2 = 4;
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);


warning('Next round'); pause
% % % % % % % % % % % % % % % % % % % Round 2 - Alpha FFC Amp versus Phase % % % % % 
% % % % % % % % % % % % % % % % % % Delay Animal L
% Animal selection
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 

sfc_mode1 = 22.4018111;
sfc_mode2 = 22.4018111;
curr_stage1 = 3;
curr_stage2 = 3;
freqband_stats1 = [10 12];
freqband_stats2 = [10 12];
plotmode1 = 1;
plotmode2 = 4;
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);

%
% % % % % % % % % % % % % % % % % % Delay Animal O
% Animal selection
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 

sfc_mode1 = 22.4018111;
sfc_mode2 = 22.4018111;
curr_stage1 = 3;
curr_stage2 = 3;
freqband_stats1 = [10 12];
freqband_stats2 = [10 12];
plotmode1 = 1;
plotmode2 = 4;
[wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude);


%% Figure 4_08 - Analyze theta peak in mode 22.4714111
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg4_08';

do_units = 0;

% SFC mode
if do_units==0
    sfc_mode = 22.4714111;
    perm_mode = 22.4714111;
%     sfc_mode = 41.6014111;
%     perm_mode = 41.6014111;
elseif do_units == 1
    sfc_mode = 52.7002010;
    perm_mode = 52.7000010;
elseif do_units == 2
    sfc_mode = 52.7000010;
    perm_mode = 52.7000010;
elseif do_units == 4        % Spectrogram
    sfc_mode = 23.4014111;
    perm_mode = 23.4014111;
end


% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 3;

if any(do_units == [1,4])
    curr_stage_sfc = 4;
end

plot_on_spect = 1;
plot_on_spectrogram = 0;
plot_on_scatter = 0;
plot_on_bargraph = 1;

switch do_units
    case 2; plot_on_spect = 0;
    case 4;
        plot_on_spectrogram = 1;
        plot_on_spect = 0;
        plot_on_bargraph = 0;
end


    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 1;
    opts_exclude.excludeL = 1;
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % % % % % Monkey L % % % % % 
    opts_exclude.excludeL = 0;
    opts_exclude.excludeO = 1; 
    
    
    perm2pls_split_plusminus=0;    % No splitting
    [wrkspc_buffer, group] = Fg_4_08_test_theta_peak(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,perm2pls_split_plusminus);
    perm2pls_split_plusminus=1;    % No splitting
    [wrkspc_buffer, group_s] = Fg_4_08_test_theta_peak(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,perm2pls_split_plusminus);
    
    inds = 1:length(group);
    if opts_exclude.excludeO == 1
        inds = [1]; inds_s = [1,2]; figure; [h1] = plot_matrix3D_custstruct([],[group(inds) group_s(inds_s)],opts_PM3D,opts_PM3Dcs);
        inds = [5]; inds_s = [9,10]; figure; [h1] = plot_matrix3D_custstruct([],[group(inds) group_s(inds_s)],opts_PM3D,opts_PM3Dcs);
    elseif opts_exclude.excludeL == 1
    end
    
    
    warning('Next round'); pause
    % % % % % Monkey O % % % % % 
    opts_exclude.excludeL = 1;
    opts_exclude.excludeO = 0; 
    
    
    perm2pls_split_plusminus=0;    % No splitting
    [wrkspc_buffer, group] = Fg_4_08_test_theta_peak(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,perm2pls_split_plusminus);
    perm2pls_split_plusminus=1;    % No splitting
    [wrkspc_buffer, group_s] = Fg_4_08_test_theta_peak(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,perm2pls_split_plusminus);
    
    inds = 1:length(group);
    if opts_exclude.excludeO == 1
        inds = [1]; inds_s = [1,2]; figure; [h1] = plot_matrix3D_custstruct([],[group(inds) group_s(inds_s)],opts_PM3D,opts_PM3Dcs);
        inds = [5]; inds_s = [9,10]; figure; [h1] = plot_matrix3D_custstruct([],[group(inds) group_s(inds_s)],opts_PM3D,opts_PM3Dcs);
    elseif opts_exclude.excludeL == 1
        inds = [2]; inds_s = [3,4]; figure; [h1] = plot_matrix3D_custstruct([],[group(inds) group_s(inds_s)],opts_PM3D,opts_PM3Dcs);
        inds = [6]; inds_s = [11,12]; figure; [h1] = plot_matrix3D_custstruct([],[group(inds) group_s(inds_s)],opts_PM3D,opts_PM3Dcs);
    end
    
    
    warning('Next round'); pause
    % With splitting
    groupmode=1;
    [wrkspc_buffer, group_s] = Fg_4_05_FFC_perm_split_plusminus(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,groupmode);
    
    if plot_on_spect
        inds = [1]; inds_s = [1,2];
        figure; [h1] = plot_matrix3D_custstruct([],[group(inds) group_s(inds_s)],opts_PM3D,opts_PM3Dcs);
        inds = [2]; inds_s = [3,4];
        figure; [h1] = plot_matrix3D_custstruct([],[group(inds) group_s(inds_s)],opts_PM3D,opts_PM3Dcs);
    end




%% Figure 4_09 - SFC phi vs ctgsetli responses
% This works, but required some fine tuning of selection parameters.
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg4_09';


% Set up pls fv (popts)
% Unit exclusion
opts_exclude = Opts_Exclude;
opts_exclude.exclude_clipping = 1;
opts_exclude.exclude_60 = 0;
opts_exclude.exclude_nans = 1;
opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC

    
sfc_mode =  2.201711101;
groupmode = 1;
sp_theshold = 20;

opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 
wrkspc_buffer = Fg_4_09_sfc_phi_vs_FFC_ctgs(wrkspc_buffer,opts_exclude, sfc_mode, groupmode,sp_theshold)

%
sfc_mode =  2.201911101;
groupmode = 1;
sp_theshold = 40;

opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 
wrkspc_buffer = Fg_4_09_sfc_phi_vs_FFC_ctgs(wrkspc_buffer,opts_exclude, sfc_mode, groupmode,sp_theshold)

%%
warning('Continuing.'); pause

sfc_mode =  2.201711101;
groupmode = 0;

opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 
wrkspc_buffer = Fg_4_09_sfc_phi_vs_FFC_ctgs(wrkspc_buffer,opts_exclude, sfc_mode, groupmode,sp_theshold)



% 

sfc_mode =  2.201911101;
groupmode = 0;

opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 
wrkspc_buffer = Fg_4_09_sfc_phi_vs_FFC_ctgs(wrkspc_buffer,opts_exclude, sfc_mode, groupmode,sp_theshold)




%% Figure 4_13f - FFC vs Unit FR for Bias (alpha / beta) with simultaneous measure of alpha and beta
% Version of Fig4_13e with expanded options. Set group.criteria values to 2
% in order to get back to Fig4_13e. Actually I just deleted Fig4_13e since
% it's redundant. You can get it by using groupmode -1 and -2.
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg4_13f';



% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 3;

    % Pls switches
    freqband_stats = [16 20];
    freqband_stats_perm = [16 20];
    freqband_stats_perm2 = [10 12];

    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 0;
    opts_exclude.excludeL = 1;
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC

        


% Animal L - Cat vs Dog
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 
groupmode=1;
do_units = 1; [wrkspc_buffer, group] = Fg_4_13f_Units_vs_FFC_simultaneous(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);
do_units = 2; [wrkspc_buffer, group] = Fg_4_13f_Units_vs_FFC_simultaneous(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);

% Animal O - Fat vs thin
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 
groupmode=2;
do_units = 1; [wrkspc_buffer, group] = Fg_4_13f_Units_vs_FFC_simultaneous(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);
do_units = 2; [wrkspc_buffer, group] = Fg_4_13f_Units_vs_FFC_simultaneous(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);


% % % % % % % % Cases 3-6; didn't work
% % % % % % Animal L - Cat vs Dog
% % opts_exclude.excludeL = 0;
% % opts_exclude.excludeO = 1; 
% % groupmode=3;
% % do_units = 1; [wrkspc_buffer, group] = Fg_4_13f_Units_vs_FFC_simultaneous(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);
% % do_units = 2; [wrkspc_buffer, group] = Fg_4_13f_Units_vs_FFC_simultaneous(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);
% % groupmode=4;
% % do_units = 1; [wrkspc_buffer, group] = Fg_4_13f_Units_vs_FFC_simultaneous(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);
% % do_units = 2; [wrkspc_buffer, group] = Fg_4_13f_Units_vs_FFC_simultaneous(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);
% % 
% % 
% % % Animal O - Fat vs thin
% % opts_exclude.excludeL = 1;
% % opts_exclude.excludeO = 0; 
% % groupmode=5;
% % do_units = 1; [wrkspc_buffer, group] = Fg_4_13f_Units_vs_FFC_simultaneous(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);
% % do_units = 2; [wrkspc_buffer, group] = Fg_4_13f_Units_vs_FFC_simultaneous(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);
% % groupmode=6;
% % do_units = 1; [wrkspc_buffer, group] = Fg_4_13f_Units_vs_FFC_simultaneous(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);
% % do_units = 2; [wrkspc_buffer, group] = Fg_4_13f_Units_vs_FFC_simultaneous(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);


%% Figure 4_13g - As 4_13f (unit biases selected by alpha/beta peaks); however, here using spectrogram perm mode
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fg4_13g';



% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 4;       % For spectrogram inputs

    % Pls switches
    freqband_stats = [16 20];
    freqband_stats_perm = [16 20];
    freqband_stats_perm2 = [10 12];

    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 0;
    opts_exclude.excludeL = 1;
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC

        


% Animal L - Cat vs Dog
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 1; 
groupmode=1;
do_units = 1; [wrkspc_buffer, group] = Fg_4_13g_Units_vs_FFC_simultaneous_spectrogram(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);
do_units = 2; [wrkspc_buffer, group] = Fg_4_13g_Units_vs_FFC_simultaneous_spectrogram(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);

% Animal O - Fat vs thin
opts_exclude.excludeL = 1;
opts_exclude.excludeO = 0; 
groupmode=2;
do_units = 1; [wrkspc_buffer, group] = Fg_4_13g_Units_vs_FFC_simultaneous_spectrogram(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);
do_units = 2; [wrkspc_buffer, group] = Fg_4_13g_Units_vs_FFC_simultaneous_spectrogram(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units);


