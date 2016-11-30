
%% Stage
% clear
% clc
run_setdefaultfig
addpath(genpath('./funcs_plots_preferred/'));

% Data loading parameters
recalc_md = 0;
recalc_sfc = 0;
recalc_ns = 0;
recalc_ev = 0;
recalc_permute = 1;
recalc_bootstrap = 0;
recalc_cfc = 0;
recalc_spect = 0;
reload_data = 0;
if get_iscromer; stage_range = [2,3];
else stage_range = [2,3];
end

% sfc_mode=2.2015100;
% perm_mode=2.201511;   % Make sure these are aligned
boot_mode = 2.201512;
    isolate_clustered_protocols = 0; % 0 for all, 1 for cells with clustered switch trials, 2 for only cells with no clustering of switch trials
spect_mode = 3.20101;
%spect_mode = 8.00001;
cfc_mode = 40.60131;

%perm_mode=sfc_mode+0.000001;   % Make sure these are aligned

gr_submode = 1;             % Which of the ctgs to plot from the bndry, RT modes (i.e. congruent/incongruent; 50 or 60% morphs

[~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
fprintf(['Chosen SFC mode: ' fname_suffix_sfc ' \n']);

% Load preferred units
    % Stage selection
    curr_stage=2;
    curr_stage_alt=3;
    curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
    curr_stage_stats = curr_stage_sfc;

    % Unit exclusion
    exclude_clipping = 1;
    exclude_60 = 0;
    exclude_lowfiring = 1;  % Justify removing low firing because of FR vs Cave plot - negative association between FR and SFC
    excludeL = 0;
    excludeO = 0;
    
    % Sensitivity switches
    sens_std_normalize = 1;     % Turning this on gives much stronger stats. It appears consistency in FR difference is more important than magnitude!!
    sort_on = 0;    % Sorting of data (this is kind of slow)
    use_morph_sens = 0;             % Loads sensitivities calculated from percent morphs, instead of using spike rate as a surrogate
    swap_in_PEV = 1;        % Use percent explained variance instead of sensitivity measures
    
    % Grouping switches
    groups_separate_monkeys = 0;
    preferred_mode = 3;     % 1 for statistics (ttest, not working in plots_preferred_unitpairs.m); 2 for raw threshold; 3 for quantile
    
    % Pls switches
    do_fisher=0;
    % freqband_stats = [40 60];
    % freqband_stats = [30 50];
    % freqband_stats = [19 39];
    freqband_stats = [20 32];
    % freqband_stats = [20 30];
    % freqband_stats = [2 4];
    % freqband_stats = [14 16];
    % freqband_stats = [18 22];
%     freqband_stats = [15 30];
    plotmode = 1;               % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8 - CFC
        ispaired = 0;
    perm2pls = 0;               % Instead of showing raw SFC values, show % cells successfully passing permutation test
        perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
        perm2pls_dophi = 0;
    
    opts_PM3Dcs.paperfig_mode=0;
    opts_PM3Dcs.do_sgolay=0;
    opts_PSC.paperfig_mode=0;
    opts_PPS.paperfig_mode = 1;
    opts_PPS.pd0 = 100;
    
    % CFC parameters
    freqband_corrcoef = [11];    % Measure correlations with this frequency
    
    

%% Sensitivity Parameters


clear senscstr
i=0;
i=i+1; senscstr{i} = 's1';
i=i+1; senscstr{i} = 's2';
% i=i+1; senscstr{i} = 's1 - sqrt(s1.*s2)';
% i=i+1; senscstr{i} = 'sqrt(s1.*s2)';
% i=i+1; senscstr{i} = '(max(s1-s2.^2,0))';
% i=i+1; senscstr{i} = '(max(s2-s1.^2,0))';
% i=i+1; senscstr{i} = 's1.*s2';
% i=i+1; senscstr{i} = 's1.^2+s2.^2';


% Regression rel vs irrel
clear sensregstr; i=0;
i=i+1; sensregstr{i} = 's1';
i=i+1; sensregstr{i} = 's2';
% i=i+1; sensregstr{i} = 's1-sqrt(s1.*s2)';
% i=i+1; sensregstr{i} = 'sqrt(s1.*s2)';
% i=i+1; sensregstr{i} = 's2-sqrt(s1.*s2)';
% i=i+1; sensregstr{i} = 'sqrt(sa1.*sa2)';
% i=i+1; sensregstr{i} = 's2.^2';
%i=i+1; sensregstr{i} = 's3';
%i=i+1; sensregstr{i} = 's4';
%i=i+1; sensregstr{i} = 's8';    % Product
% i=i+1; sensregstr{i} = 'sa1';
% i=i+1; sensregstr{i} = 'sa2';
%i=i+1; sensregstr{i} = 'sa3';
%i=i+1; sensregstr{i} = 'sa4';
%i=i+1; sensregstr{i} = 'sa8';    % Product
% i=i+1; sensregstr{i} = 's9';    % Firing rate
i=i+1; sensregstr{i} = 'sa9';    % Alternate firing rate (for delay stage, most likely - we should use this)
% i=i+1; sensregstr{i} = 'cn5';    % Percent switch trials
% i=i+1; sensregstr{i} = 'cn1';       % Monkey dates
% i=i+1; sensregstr{i} = 'cn2';       % Monkey L dates
% i=i+1; sensregstr{i} = 'cn3';     % Monkey O dates
% i=i+1; sensregstr{i} = 'cn4';     % Monkey L or O (binary)
% i=i+1; sensregstr{i} = 'spc1.*cn1';
% i=i+1; sensregstr{i} = 'spc2.*cn2';
% i=i+1; sensregstr{i} = 's1.*cn4';
% i=i+1; sensregstr{i} = 's2.*cn4';
% i=i+1; sensregstr{i} = 'sa1.*cn4';
% i=i+1; sensregstr{i} = 'sa2.*cn4';
% i=i+1; sensregstr{i} = 's1.*~cn4';
% i=i+1; sensregstr{i} = 's2.*~cn4';
% i=i+1; sensregstr{i} = 'sa1.*~cn4';
% i=i+1; sensregstr{i} = 'sa2.*~cn4';

if (sum(strcmp(sensregstr,'cn2')) > 0 || sum(strcmp(sensregstr,'cn3')) > 0 ) && sum(strcmp(sensregstr,'cn4')) == 0
    warning('cn2 or cn3 included in regression, but cn4 is not. Cn4 must be including in order to make cn2 and cn3 scale invariant.');
        % This is an issue because cn2 and cn3 are forced to be zero zeros for monkeys O and
        % L, respectively. Therefore need to introduce an additional
        % constant shift to allow these to vary truly independently.
    keyboard;
end

% Dependent variable
i=i+1; sensregstr{i} = 'p11';

% % Regression sample vs delay
% clear sensregstr; i=0;
% i=i+1; sensregstr{i} = 's8';
% i=i+1; sensregstr{i} = 'sa8';
% i=i+1; sensregstr{i} = 's9';

sensregnames = sensregstr;
sensregnames = strrep(sensregnames,'.*','');
sensregnames = strrep(sensregnames,'sqrt(','');
sensregnames = strrep(sensregnames,')','');

% sensregnames = {'SampleA','SampleB','FRate','Date','SFC'}




%% Plot Switches

% General
domean=1;
stats_dozscore = 0;
do_zscore = 0;
showError = 1;

% Plotting toggles
% Time series
    plot_SFC_spectra = 0;
    
% Stats
    plot_SFC_spectra_stats = 0;
    fig_plot_histplot = 0;
    fig_plot_venn = 0;

% Spectrograms
    plot_spectrogram = 0;
        caxismax = 200;
        plot_spectrogram_vs_ev = 0;
        plot_spectrogram_individuals = 0;       % Spectrograms for each individual cell
        spectrogram_thresh = 0;                 % Plots percentage of cells above a given SFC threshold
        
% Evokedes
    plot_evoked = 0;
        plot_ev_raw = 1;
        plot_ev_psd = 1;
        plot_ev_spect = 0;

    plot_badfiles_venn = 0;

% Scatter plots
    scattergroups_sens_vs_sens = 0;
    scattergroups_sens_vs_pls = 0;
    regress_pls_vs_sens = 0;
    plot_fisher_test_suite = 0;
    
    
% PCA
    plot_PCA = 0;
    
% Surverys
    find_archetypal_beta_cells = 0;
    plot_find_high_beta_cells = 0;
    phi_survey_ctgs = 0;
    survey_boundary_vs_categorization_index = 0;     % For boundary trials
    Cave_survey_ctgs = 0;        % Like Phi survey, searches for significantly different Ca across ctgs
    
% Permutation test results
    plot_Sch_sens = 0;   % Does the same thing as above, but with more specialized code
    bootstrap_validation_KStest = 0;        % Makes sure data from bootstrap mode is semi-Gaussian, for bootstrapt est

% For dealing with phase plots
    if plotmode == 4
        %ispaired = 1;
    end
    

    
%% Groups Definition
% Setup groups

symmetry_mode = 3;          % 1-A; 2-B; 3-By Preference; 4-TestSymmetry_Preferred; 5-TestSymmetry_NonPreferred

% Swap data for test stage, if that's where we're looking
% Activate swap mode
swap_mode = 0;
use_alternate_stage = 0;

switch swap_mode
    case 0;
        swap_map = [1,2,3,4 5 6 7 8; [1,2,3,4 5 6 7 8] ];   % Identity (only if we're using swap and not move); THIS BUGGERS UP WHEN USING CTG 11
    %case 1; swap_map = [1,2,3,4; [9,10,11,12] ];    % Test sensitivity; maps ctgs 1,2,3,4 to 9,10,11,12
    %case 2; swap_map = [1,2,3,4; [13,14,15,16] ];   % match non-match, schA schB
    %case 3; swap_map = [1 2; [13 14] ];             % match non-match
    case 4; swap_map = [1 2; [9 10] ];             % schA schB
end
swap_map_pref = unique(round(swap_map/2)','rows')'; % Maps the ctg indices onto pref indicies (i.e. src -> spc). For example 1-4 maps onto 1-2; etc.

% [gr,mask] = group_setup_FR4D;                      % Generate grouping structure
[gr,mask] = group_setup_FR2D;       % Generate grouping structure
[gr.all_cells, gr.all_sens, gr.all_sens4, gr.single, gr.single_sens] = group_setup_other;

% Scored by Right minus wRong
%gr.all = normalize_diffs(gr.all,'diff');
gr.rmr = score2grouping(gr.all,[gr.all.diff]);                     % Don't use Mask for rmr
% gr.rmr = score2grouping(gr.all(mask),[gr.all(mask).diff]);         % Use Mask for rmr
% gr_netpref = score2grouping(gr.all,[gr.all.netprefA]); gr_netpref.ctgs = []

gr = group_enumerate_sens(gr,mask,symmetry_mode);     % Enumerate grouping structure across various "sensitivity tests"
gr = group_enumerate_stages(gr);                      % Create a grouping structure that can handle compairsons between two stages
[gr.SD_QP, gr.S_QP] = group_enumerate_all(gr.all,gr.single);     % Create a group that compares all stages and all inputs (P-preferred and Q-Cue)
gr.sd_newcoords = group_coords2to4(gr.sd);
group = gr.all;

% group = gr.all([1:4 11]);
% group = gr.all([1 5 9 13 11]);
% group = gr.all([1 2 5 6]);                             % Q1
% group = gr.all([3 4 7 8]); group = gr.all([3 4]);      % Q2
% group = gr.all([9 10 13 14]); group = gr.all([9 13]);  % Q3
% group = gr.all([11 12 15 16]); group = gr.all([11 12 15]); % Q4
% group(end+1) = gr.all(11);

% %% Test cue classifiers vs others
% clear group;
% i=1; group(i) = group_merge(gr.SD_QP(1:8)); group(i).legend = 'P1P2';
% i=2; group(i) = group_merge(gr.SD_QP(13:end-1)); group(i).legend = 'Q1';
% i=3; group(i) = group_merge(gr.SD_QP(end));group(i).legend = 'Q0P0';

% %% Sample vs delay neurons P1
% clear group
% i=1; group(i) = group_merge(gr.sd_Dall(4:6)); group(i).legend = 'Sample P1';
% i=2; group(i) = group_merge(gr.sd_Dall([2,5,8])); group(i).legend = 'Delay P1';

% %% Sample vs delay neurons P2
% clear group
% i=1; group(i) = group_merge(gr.sd_Dall(1:3)); group(i).legend = 'Sample P2';
% i=2; group(i) = group_merge(gr.sd_Dall([1,4,7])); group(i).legend = 'Delay P2';

% %% Multitaskers
% clear group
% i=1; group(i) = group_merge(gr.all([1,3])); group(i).legend = 'SchA';
% i=2; group(i) = group_merge(gr.all([1,2])); group(i).legend = 'SchB';


% For Venn
% clear group
% group(1:2) = Grp;
% group(1) = group_merge(gr.all(1),gr.all(3)); group(1).legend = ('Sch A');
% group(2) = group_merge(gr.all(1),gr.all(2)); group(2).legend = ('Sch B');

% For FR traces
% group(1).legend = ('Sch A Preferred');
% group(2).legend = ('Sch A Non-Preferred');
%
% Rearrange for gr.all
% group(2:3) = [group(3) group(2)];

% group(1).legend = ('Multitaskers');
% group(2).legend = ('Specialists');
% group(3).legend = ('No Preference');
% 
% For Sens_Sch
% group(1).legend = 'Specialist Neurons, Preferred Sch';
% group(2).legend = 'Specialist Neurons, Non-Preferred Sch';


if swap_mode ~= 0
    group = group_swap.crit (group,swap_map_pref);      % THIS BUGGERS UP WHEN USING CTG 11
    group = group_swap.ctgs (group,swap_map);
end
if use_alternate_stage; group = group_swap.alt (group); end


if get_iscromer; group = group_swap.cromer_shift(group); end

if groups_separate_monkeys
    group_M(1:2) = Grp;
    group_M(1).criteria_sfc(2) = 1; group_M(1).legend = 'Monkey1';
    group_M(2).criteria_sfc(2) = 0; group_M(2).legend = 'Monkey2';
    
    %group = group_OUTER(group,group_M); % Separate all groups according to Monkeys
    group = group_M;        % Plot just 2 separate Monkeys
end



N_criteria = 5;

% Override group special modes (boundary trials, reaction times, switch trials)
group_override = get_group_overrides(fname_suffix_sfc, N_criteria,gr_submode);
if ~isempty(group_override); group = group_override; end


group.getcrit;
% group.gettable('ctgs');

%% Import data if running as a function; define dependent parameters

running_as_script = ~is_calledby;
if ~running_as_script
    pop_workspace(false);
end

[~, perm_subgroups] = decode_sfc_mode(perm_mode); [fname_suffix_perm] = build_sfcmode(perm_mode, perm_subgroups);

    
%% Load data    
% Generate some basic parameters and populate groups structure
fs = 1/get_dt;

% Estimate file list and file_ranges
% [file_list, file_range] = estimate_filelist(sfc_mode, 2);
file_range = 1:79;

path_matfiles_store = getpath('path_matfiles_store');


% Create buffer variable
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end

% Load Metadata
buff_fieldname = 'currmd';
buff_func = @() func_recalc_md(file_range);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_md,reload_data);
md = wrkspc_buffer.(buff_fieldname).md;


% Load the 2D data (coherence)
% SFC
buff_fieldname = ['sfc_' mode2modename(sfc_mode)];
buff_func = @() func_recalc_sfc(stage_range,file_range,sfc_mode);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_sfc,reload_data);
vars_pull(wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc)));

% Percent trials switch (used for sens, but it's really not a useful predictor)
buff_fieldname = ['switchtr'];
buff_func = @() calc_percent_switch(md);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_sfc,reload_data);
vars_pull(wrkspc_buffer.(buff_fieldname));

% Load specialization indices (statsFR)
load(fullfile(getpath('path_buffer_specialization_stats'),'specialization.mat'));
statsnames = stats.(stagename(curr_stage_stats)).stats(1).m.ctgsetnames;
statsFR = stats.(stagename(curr_stage_stats)).mu_arr;
statsFR = squeeze(statsFR);

% Load permutation tests
% buff_mode = 22.201311;
% buff_mode = 22.401311;
buff_mode = perm_mode;
buff_fieldname = ['per_' mode2modename(buff_mode)];
buff_func = @() func_recalc_perms(buff_mode,stage_range,file_range);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_permute,reload_data);
% vars_pull(wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc)));

% Load bootstrap tests
buff_mode = boot_mode;
buff_fieldname = ['boot_' mode2modename(buff_mode)];
buff_func = @() func_recalc_boot(buff_mode,[3],file_range);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_bootstrap,reload_data);
% vars_pull(wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc)));



% Load ns (spike trains)
buff_mode = 5.1;
buff_fieldname = ['ns_' mode2modename(buff_mode)];
buff_func = @() load_and_merge(file_range, buff_mode, 5, {'spikerates'});
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_ns,reload_data);
vars_pull(wrkspc_buffer.(buff_fieldname));
ns = spikerates; clear spikerates

% Load evoked (mean LFP)
buff_mode = 9.00001;
buff_fieldname = ['ev_' mode2modename(buff_mode)];
buff_func = @() load_and_merge(file_range, buff_mode, 5, {'E'});
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_ev,reload_data);
vars_pull(wrkspc_buffer.(buff_fieldname));
ev = E; clear E;    % Rename E to ev

% Load CFC data
if curr_stage_sfc == 3
    buff_fieldname = ['cfc_' mode2modename(cfc_mode)];
    buff_func = @() func_recalc_mode40 (cfc_mode,stage_range,file_range);
    wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_cfc,reload_data);
    %vars_pull(wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc)));
    s = wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc));
    Ptr = s.Ptr;
    Rtr = s.Cave;
    ftr = s.f(1,:);
    clear s
end


% Load 3D data (coherograms)
buff_fieldname = ['ssp_' mode2modename(spect_mode)];
buff_func = @() func_recalc_3D (spect_mode,file_range);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_spect,reload_data);
vars_pull(wrkspc_buffer.(buff_fieldname));
T = T-2;

% Setup some metadata variables
N_ctgs_base = N_ctgs_base(1);
N_ctgs_extras = N_ctgs_extras(1);
N_ctgs = N_ctgs_base + N_ctgs_extras+1;

% Bad files
[bad_lowfiring] = get_lowfiring(spikerate_mu);
[bad_clip, bad_60,fnames_clip] = loadall_bad_data(curr_stage_badclipping,file_range,md);
