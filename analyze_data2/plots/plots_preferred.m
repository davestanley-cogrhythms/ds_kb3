
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
recalc_permute = 0;
recalc_bootstrap = 0;
recalc_cfc = 0;
recalc_spect = 0;
reload_data = 0;
if get_iscromer; stage_range = [2,3];
else stage_range = [2,3];
end

sfc_mode=2.2015100;
perm_mode=2.2015111;   % Make sure these are aligned
boot_mode = 2.201512;
    isolate_clustered_protocols = 0; % 0 for all, 1 for cells with clustered switch trials, 2 for only cells with no clustering of switch trials
spect_mode = 3.20101;
%spect_mode = 8.00001;
cfc_mode = 40.60131;

%perm_mode=sfc_mode+0.000001;   % Make sure these are aligned

gr_submode = 1;             % Which of the ctgs to plot from the bndry, RT modes (i.e. congruent/incongruent; 50 or 60% morphs


% Load preferred units
    % Stage selection
    curr_stage=3;
    curr_stage_alt=3;
    curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
    curr_stage_stats = curr_stage_sfc;

    % Unit exclusion
    exclude_clipping = 1;
    exclude_60 = 0;
    exclude_lowfiring = 0;  % Justify removing low firing because of FR vs Cave plot - negative association between FR and SFC
    excludeL = 0;
    excludeO = 0;
    
    % Sensitivity switches
    sens_std_normalize = 1;     % Turning this on gives much stronger stats. It appears consistency in FR difference is more important than magnitude!!
    sort_on = 0;    % Sorting of data (this is kind of slow)
    use_morph_sens = 0;             % Loads sensitivities calculated from percent morphs, instead of using spike rate as a surrogate
    swap_in_PEV = 1;        % Use percent explained variance instead of sensitivity measures
        use_sqrt_of_biased_PEV = 1;
    
    % Grouping switches
    groups_separate_monkeys = 0;
    preferred_mode = 1;     % 1 for statistics (ttest, not working in plots_preferred_unitpairs.m); 2 for raw threshold; 3 for quantile
    
    % Pls switches
    do_fisher=0;
    % freqband_stats = [40 60];
    % freqband_stats = [30 50];
%     freqband_stats = [19 39];
%     freqband_stats = [20 32];
    freqband_stats = [16 20];
    % freqband_stats = [2 4];
    % freqband_stats = [14 16];
    % freqband_stats = [18 22];
%     freqband_stats = [15 30];
    plotmode = 1;               % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8 - CFC
        ispaired = 0;
    permdat2pls = 0;            % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    perm2pls = 0;               % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
        perm2pls_dophi = 0;
        sort_pls = 0;           % Sort into preferred and non-preferred
        do_diff = 0;            % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
    
    opts_PM3Dcs.paperfig_mode=0;
    opts_PM3Dcs.do_sgolay=0;
    opts_PSC.paperfig_mode=0;
    opts_PPS.paperfig_mode = 1;
    opts_PPS.pd0 = 100;
    
    % CFC parameters
    freqband_corrcoef = [11];    % Measure correlations with this frequency
    

[~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
fprintf(['Chosen SFC mode: ' fname_suffix_sfc ' \n']);

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
i=i+1; sensregstr{i} = 'p10';

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
    plot_SFC_spectra = 1;
    
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
    scattergroups_sens_vs_sens = 1;
    scattergroups_sens_vs_pls = 0;
    regress_pls_vs_sens = 1;
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
    
%% Import data if running as a function; define dependent parameters

running_as_script = ~is_calledby;
if ~running_as_script
    pop_workspace(false);
end

[~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
[~, perm_subgroups] = decode_sfc_mode(perm_mode); [fname_suffix_perm] = build_sfcmode(perm_mode, perm_subgroups);

    
%% Groups Definition
if ~exist('group','var')
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
    gr2 = gr.all_cells;             % Alias
    % Scored by Right minus wRong
    %gr.all = normalize_diffs(gr.all,'diff');
    gr.rmr = score2grouping(gr.all,[gr.all.diff]);                     % Don't use Mask for rmr
    % gr.rmr = score2grouping(gr.all(mask),[gr.all(mask).diff]);         % Use Mask for rmr
    % gr_netpref = score2grouping(gr.all,[gr.all.netprefA]); gr_netpref.ctgs = []

    gr = group_enumerate_sens(gr,mask,symmetry_mode);     % Enumerate grouping structure across various "sensitivity tests"
    gr = group_enumerate_stages(gr);                      % Create a grouping structure that can handle compairsons between two stages
    [gr.SD_QP, gr.S_QP] = group_enumerate_all(gr.all,gr.single);     % Create a group that compares all stages and all inputs (P-preferred and Q-Cue)
    gr.sd_newcoords = group_coords2to4(gr.sd);
    group = gr.rmr;
%     group(1:3) = group([3,2,1]);

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
    if permdat2pls || perm2pls
        group_override = get_group_overrides2(fname_suffix_sfc, N_criteria,gr_submode);
        if ~isempty(group_override); group = group_override; end
    end


    group.getcrit;
    % group.gettable('ctgs');
end
    
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

% % Load specialization indices (statsFR)
% load(fullfile(getpath('path_buffer_specialization_stats'),'specialization.mat'));
% statsnames = stats.(stagename(curr_stage_stats)).stats(1).m.ctgsetnames;
% statsFR = stats.(stagename(curr_stage_stats)).mu_arr;
% statsFR = squeeze(statsFR);

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
    Rtr = s.Cave;5
    ftr = s.f(1,:);
    clear s
end


% Load 3D data (coherograms)
% buff_fieldname = ['ssp_' mode2modename(spect_mode)];
% buff_func = @() func_recalc_3D (spect_mode,file_range);
% wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_spect,reload_data);
% vars_pull(wrkspc_buffer.(buff_fieldname));
% T = T-2;

% Setup some metadata variables
N_ctgs_base = N_ctgs_base(1);
N_ctgs_extras = N_ctgs_extras(1);
N_ctgs = N_ctgs_base + N_ctgs_extras+1;

% Bad files
[bad_lowfiring] = get_lowfiring(spikerate_mu);
[bad_clip, bad_60,fnames_clip] = loadall_bad_data(curr_stage_badclipping,file_range,md);

%% PLS variable setup
if std(phiave) > 2*pi   % Convert phiave to radians if not already. Std should be less than 2 pi. If not, convert
            phiave = phiave/180*pi; end
switch plotmode
    case {1}
        pls = Cave;
        if do_fisher; pls = fisher_trans(pls); end
        group(1).data_name =  'SFC';
        if exist('Ssp','var'); SPG = Ssp; else; fprintf('Spectrogram data not loaded \n'); end
        abscissa = f;
        
        % Stupid clause to expand pls with fake Infs if it's shorter than desired size (Nctgs = 11)
        if perm2pls || permdat2pls
            pls = pls_pad(pls,sfc_subgroups,1); % Pad 11 with 0s
        else
            pls = pls_pad(pls,sfc_subgroups,2); % Pad 4 --> 11
        end
    case 2
        pls = S1ave;
        group(1).data_name =  'LFP PSD';
        abscissa = f;
    case 3
        pls = ns;
        group(1).data_name =  'SPK s-1';
        t = [1:size(ns,1)]*get_dt;
        t = t - 2;  % Centre to sample on time
        abscissa = t;
    case 4
        pls = Cave.*exp(1i*phiave);    % Imaginary vector containing amp and phase
        pls = pls_complex2angle_centpi(pls);
        group(1).data_name =  'Phase';
        abscissa = f;
    case 5
        pls = S2ave;
        group(1).data_name =  'SPK PSD';
        SPG = S2AVE;
        abscissa = f;
    case 6
        pls = Cave.*exp(1i*phiave);    % Imaginary vector containing amp and phase
        group(1).data_name =  'ComplexMag';
        fcomplex2real = @abs;
        abscissa = f;
    case 7
        pls = Cave.*exp(1i*phiave);    % Imaginary vector containing amp and phase
        fcomplex2real = @pls_complex2angle_centpi;
        group(1).data_name =  'ComplexPhase';
        abscissa = f;
    case 8
        pls = Cave.*exp(1i*phiave);    % Imaginary vector containing amp and phase
        group(1).data_name =  'Complex';
    case 9
        ind = find(ftr >= freqband_corrcoef,1,'first');
        pls = squeeze(Rtr(ind,:,:,:));
        ulpairs = build_unit_LFP_pairs_all(md)';
        pls = pls(:,ulpairs(:,2),:);
        abscissa = ftr;
        group(1).data_name =  'CFC';
        
        % if pls is small, extend it to the full 11 dimensions
        if size(pls,3) == 2
            temp = ones(size(pls(:,:,end))) * Inf;
            pls = cat(3,repmat(temp,[1,1,10]),pls(:,:,end));
        end
        
end



% Set up axis limits
switch plotmode
    case {1}
        for i=1:length(group); if isempty(group(i).xlims_desired); group(i).xlims_desired = [0 160]; end; end
        if sfc_mode >= 2.0 && sfc_mode <= 2.2;
            if isempty(group(1).ylims_desired); group(1).ylims_desired = [0.02 0.14]; end
        end
    case {2,4:9}
        for i=1:length(group); if isempty(group(i).xlims_desired); group(i).xlims_desired = [0 90]; end; end
end




% Pls stats - we need to recalculate this later after sorting!


% ID bad data
if plot_spectrogram || plot_spectrogram_vs_ev
    bad_sfc_nan = get_bad_SFC(SPG); %Spectrogram
    curr_stage_badclipping = 3;     % Use delay stage for bad files, since this one had highest amplitude and most clipping... could alternatively use stage 5, but this is definitely over-conservative (try running: sum([bad_clip' bad_60']))
elseif plot_evoked
    bad_sfc_nan = get_bad_SFC(ev);
else
    bad_sfc_nan = get_bad_SFC(pls);
end

% Swap in switch trials if using them
[bad_switch] = get_bad_switch(sum_ctgsetli, sfc_subgroups, isolate_clustered_protocols);



%% Remove bad cells & Sort

% Get bad files
bad_any = false(size(bad_clip));
bad_any = bad_any | bad_switch;
bad_any = bad_any | bad_sfc_nan;
if exclude_clipping; bad_any = bad_any | bad_clip; end
if exclude_60; bad_any = bad_any | bad_60; end
if exclude_lowfiring; bad_any = bad_any | bad_lowfiring; end

ismonkeyL_temp=~cellfun(@isempty,(strfind(funames1D,'L')));
if excludeL; bad_any = bad_any | ismonkeyL_temp; end
if excludeO; bad_any = bad_any | (~ismonkeyL_temp); end
clear ismonkeyL_temp
    
if plot_badfiles_venn; figl; plot_venn_badfiles(bad_lowfiring, bad_clip, bad_60, bad_any); end

%% Swap in perm_mode Cave1 and Cave2 to pls
if permdat2pls
    pls = permCave2pls_swap(wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
        perm2pls_dophi);
    sfc_subgroups = perm_subgroups; fname_suffix_sfc = fname_suffix_perm;
        
    group = group_enum_ctgs(gr2,size(pls),fname_suffix_perm,[0 80]);

end


%% Swap in permutation results to pls
if perm2pls         % Instead of showing raw SFC values, show % successfully passing permutation test
    
    %Differences normalized by null distribution
    pls = perm2pls_swap(wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
        perm2pls_dophi,perm2pls_do_bh,bad_any);
    sfc_subgroups = perm_subgroups; fname_suffix_sfc = fname_suffix_perm;

    %Rework group definition (cheesy)
    group = group_enum_ctgs_diff(gr2,size(pls),fname_suffix_perm,[0 200]);
    
end


%% Sort pls into preferred and non-preferred
if sort_pls
    for i = 1:floor(size(pls,3)/2)
        first = i*2-1;
        second = i*2;
        [pls(:,:,first),pls(:,:,second)] = sort_pref_nonpref(f,pls(:,:,first),pls(:,:,second),freqband_stats);
    end
end

%% Difference between adjacent PLS
if do_diff  % pls(:,:,2)-pls(:,:,1), etc.
    pls = -1*diff(pls,[],3);
    pls = pls(:,:,1:2:end);
    if perm2pls_dophi; pls = get_circ_diff(pls); end
    
    group = group_enum_ctgs_diff(gr2,size(pls),fname_suffix_sfc, [0 80]);
end

%% Load preferred schemes, spiking rates; remove bad files
[spc,src,ssc,pr_fnames,pr_funames1D] = extract_spsr(curr_stage, bad_any);
[spa,sra,ssa,pr_fnames,pr_funames1D] = extract_spsr(curr_stage_alt, bad_any);

% Calculate percent explained variation
[pevc,eta2c] = extract_PEV(curr_stage);
[peva,eta2a] = extract_PEV(curr_stage_alt);

if use_sqrt_of_biased_PEV
    warning ('Using biased PEV');
    pevc = eta2c;
    peva = eta2a;
    % % warning('Shifting pevc and peva to be > 0 for numerical reasons; this may be incorrect');
    % % pevc = pevc - min(pevc(:));
    % % peva = peva - min(peva(:));
    % % pevc(pevc<0) = 0;
    % % peva(peva<0) = 0;
    pevc=sqrt(pevc);
    peva=sqrt(peva);
end

% Check that the correct files were loaded in both cases
    if sum(~strcmp(pr_funames1D,funames1D)) ~= 0; fprintf(['File mismatch! Likely some files are missing for either '...
        'spikerates or sfc. Aborting \n']); end
    if sum(~strcmp(pr_fnames,fnames_clip)) ~= 0; fprintf(['File mismatch! Likely some files are missing for either '...
        'spikerates or sfc. Aborting \n']); end


% Sort data so that the high firing rate categories always come first. Do this based on the CURRENT STAGE, not the alternate
if sort_on
    tic
    [src, ssc, pls, ns] = sortall_data(src,ssc,pls,ns,N_ctgs_base,N_ctgs_extras);
    toc
end



%% Sensitivities
% Load sensitivity estimates based on morphs or firing rates

[sens, sensnames] = calc_sens(src,ssc,curr_stage,use_morph_sens,sens_std_normalize);
[sens_alt, sensnames_alt] = calc_sens(sra,ssa,curr_stage_alt,use_morph_sens,sens_std_normalize); sensnames_alt = strcat('A', sensnames_alt);

if swap_in_PEV
    sens(:,1:5) = pevc;
    sens_alt(:,1:5) = peva;
%     sens(:,1:5) = eta2c;
%     sens_alt(:,1:5) = eta2a;
end

% % Cell numbers (not used)
% MLcellnumi=~cellfun(@isempty,(strfind(funames1D,'L')));MLcellnum=double(MLcellnumi);
% MLcellnum(MLcellnumi)=cumsum(MLcellnumi(MLcellnumi));
% 
% MOcellnumi=~cellfun(@isempty,(strfind(funames1D,'O'))); MOcellnum=double(MOcellnumi);
% MOcellnum(MOcellnumi)=cumsum(MOcellnumi(MOcellnumi));
% cellnum = [[MLcellnum(:) + MOcellnum(:)] MLcellnum(:) MOcellnum(:) MLcellnumi(:) percent_switch(:)];

% Cell dates
if ~get_iscromer
    ismonkeyL=~cellfun(@isempty,(strfind(funames1D,'L')));
else
    ismonkeyL=~cellfun(@isempty,(strfind(funames1D,'li'))); % Or something like that...
end
date = cell2mat({md.datenum_numeric});
Nunits = arrayfun(@(s) length(s.unit_names),md);
date = replicate_inline(date,Nunits); % Make date be 1x549 instead of 1x79 (one entry for each unit)
date = date-min(date);

dateL = date; dateL = dateL - min(date(ismonkeyL)) + 0; dateL(~ismonkeyL) = 0;
dateO = date; dateO = dateO - min(date(~ismonkeyL)) + 0; dateO(ismonkeyL) = 0;
cellnum = [date(:) dateL(:) dateO(:) ismonkeyL(:) percent_switch(:)];

%% Derive new spc and spa
switch preferred_mode
    case 1
        % Uses the default spc, extracted above from run_preferred3.m
    case 2
        % Original threshold method
        thresh = 0.15;
        sensrange=1:5;
        spcd = sens(:,sensrange) > thresh;
        spad = sens_alt(:,sensrange) > thresh;
        spc=spcd;
        spa=spad;
    
    case 3
        % Quantile thresholding method
        quant_thresh = 0.7;
        sensrange=[1:4, 10];

        N = size(sens,1);
        sens_curr = sens; temp = quantile(sens_curr(~bad_any,sensrange),quant_thresh); spcd = [sens_curr(:,sensrange) > repmat(temp,N,1)];
        sens_curr = sens_alt; temp = quantile(sens_curr(~bad_any,sensrange),quant_thresh); spad = [sens_curr(:,sensrange) > repmat(temp,N,1)];
        clear temp
        spc = spcd;
        spa = spad;
    
end


%% Parse sensitivity strings used for scatterplots and regression
if sfc_subgroups(2) == 0 && (scattergroups_sens_vs_sens | scattergroups_sens_vs_pls | regress_pls_vs_sens)
    pls_stats_temp = calc_pls_stats(abscissa,pls,freqband_stats,'do_mean_ctgs',0,'do_mean_freq',0);
    [sens_sc bad_zeros_sc] = parse_sens(senscstr,sens,sens_alt,spc,spa,pls_stats_temp,cellnum);
    [sens_reg bad_zeros_reg]= parse_sens(sensregstr,sens,sens_alt,spc,spa,pls_stats_temp,cellnum);
    clear pls_stats_temp
end



%% Quartile calculations
% Quartile: Splits the last group into SFC and non-SFC cells
do_quartile = 0;

quartile_mode = 1;
switch quartile_mode

    case 1;                                                                 % Based on Q4
        pls_stats_temp = calc_pls_stats(f,Cave(:,:,end),freqband_stats,'do_mean_ctgs',1);

        q4 = quantile(pls_stats_temp,0.75);
        q4_sfc = pls_stats_temp >= q4;
        clear pls_stats_temp
        legendtext1 = 'SFC'; legendtext2 = 'Non-SFC';

    case 2;
        q4_sfc = wrkspc_buffer.sfc_15_2.stage3.p(1,:,end) <= 0.05;       % Based on glm statistics
        legendtext1 = 'SFC'; legendtext2 = 'Non-SFC';
        sum(q4_sfc)
    case 3                                                                  % Based on Chronux mag coherence
        s = wrkspc_buffer.sfc_2_2011.stage3;                        
        q4_sfc = calc_significant_chronux(s.f,s.Cave,s.confC,mean(freqband_stats));
        q4_sfc = q4_sfc(:,end); 
%         q4_sfc = q4_sfc(:,9) & q4_sfc(:,10);
        sum(q4_sfc)
        legendtext1 = 'SFC'; legendtext2 = 'Non-SFC';
    case 4                                                                  % Based on Chronux Jackknife
        s = wrkspc_buffer.sfc_2_2012.stage3;                        
        q4_sfc = calc_significant_chronux_jk(s.f,s.Cave,s.Cerr1,mean(freqband_stats));
        %q4_sfc = q4_sfc(:,end); 
        q4_sfc = q4_sfc(:,9) & q4_sfc(:,10);
        sum(q4_sfc)
        legendtext1 = 'SFC'; legendtext2 = 'Non-SFC';
    case 5
        s = wrkspc_buffer.sfc_2_2012.stage3;                        % Based on Chronux phi difference
        %[q4_sfc1, z1, zstd1, t1] = calc_significant_chronux_phi_bh_stepup(s.f,s.phiave,s.phistd,0.05,mean(freqband_stats),[1 2],1);
        %[q4_sfc2, z2, zstd2, t2] = calc_significant_chronux_phi_bh_stepup(s.f,s.phiave,s.phistd,0.05,mean(freqband_stats),[3 4],1); 
        %[q4_sfc3, z3, zstd3, t3] = calc_significant_chronux_phi_bh_stepup(s.f,s.phiave,s.phistd,0.05,mean(freqband_stats),[9 10],1);
        do_bh=1;
        [q4_sfc , z , zstd , t ] = calc_significant_chronux_phi_bh_stepup(s.f,s.phiave,s.phistd,0.05,bad_any,mean(freqband_stats),[9 10],do_bh); sum(q4_sfc)
        legendtext1 = 'Phi Sch Sensitive'; legendtext2 = 'Phi Sch Insensitive';
    case 6
        s = wrkspc_buffer.sfc_2_2012.stage3;   
        do_bh=0;
        [q4_sfc ,pvals] = calc_significant_chronux_delCa(s.f,s.Cave,s.Cerr1,s.Cerr2,0.05,bad_any,mean(freqband_stats),[9 10],do_bh); sum(q4_sfc)
        legendtext1 = 'SFC Sch Sensitive'; legendtext2 = 'SFC Sch Insensitive';
    case 7                              % Cave permutation test for SchA vs SchB
        s = wrkspc_buffer.per_2_201511.stage3;
        do_bh = 1;
        %[q4_sfc, pvals1] = calc_significant_permutation_old_withplots(s.f,s.dCave,s.Cavep1,s.Cavep2,0.05,bad_any,mean(freqband_stats),do_bh,Cave); sum(q4_sfc)
        [q4_sfc, pvals1] = calc_significant_permutation(s.f,s.pvals_Cave,0.05,bad_any,mean(freqband_stats),do_bh); sum(q4_sfc)
        legendtext1 = 'SFC Sch Sensitive'; legendtext2 = 'SFC Sch Insensitive';
    case 8                              % Phi permutation test for SchA vs SchB
        s = wrkspc_buffer.per_2_201111.stage3;
        do_bh = 1;
        %[q4_sfc, pvals1] = calc_significant_permutation_old_withplots(s.f,s.dCave,s.Cavep1,s.Cavep2,0.05,bad_any,mean(freqband_stats),do_bh,Cave); sum(q4_sfc)
        [q4_sfc, pvals1] = calc_significant_permutation(s.f,s.pvals_phi,0.05,bad_any,mean(freqband_stats),do_bh); sum(q4_sfc)
        legendtext1 = 'SFC Sch Sensitive'; legendtext2 = 'SFC Sch Insensitive';
end
% figure; plot(zscore(pvals))
% hold on; plot(zscore(pvals1),'r')

%sum(q4_sfc(:) & ~bad_any(:))

%bad_any = bad_any | ~q4_sfc;

q4_sfc = q4_sfc(:);
if do_quartile
    Ngr = length(group);
    group(Ngr+1) = group(Ngr);
    i=Ngr;group(i).criteria_sfc = repmat([1 2],[size(group(i).criteria,1),1]); group(i).legend=([group(i).legend ' ' legendtext1]);
    i=Ngr+1;group(i).criteria_sfc = repmat([0 2],[size(group(i).criteria,1),1]); group(i).legend=([group(i).legend ' ' legendtext2]);
    Ngr = length(group);
end


do_multitaskers_and_nondeciders = 0;
if do_multitaskers_and_nondeciders
    %group(2:end+1) = group(1:end);
    %group(1) = gr.rmr(1);
    group(end+1) = gr.rmr(3);
end

%% Calculate group cells
% Extract units sensitive to given schemes
for i = 1:length(group)
    [group(i).cells, group(i).numcells]= get_grouped_cells(group(i),[spc,spa,q4_sfc, ismonkeyL']);
end

% Remove bad files
spc(bad_any,:) = false;
spa(bad_any,:) = false;
[group] = remove_bads_cells(bad_any,group);                 % Remove bad cells from cell lists
%[bad_any,src,sra,pls,ns,ev,SPG] = remove_bads_data(bad_any,src,sra,pls,ns,ev,SPG);  % Fill dataset values with nans for cells we know are bad!
% group.gettable('numcells');
% group.gettable('right');


%% Groups Loading Setup
% Generate SFC data for plotting and statistics
for i = 1:length(group)
    switch plotmode
        case {1,2,4,5,6,7,8}
            [group(i).data, group(i).funames] = get_grouped_data(group(i).ctgs,group(i).cells,pls,funames1D);
            [group(i).datastats, group(i).freqband_stats] = calc_pls_stats(abscissa,group(i).data,freqband_stats,'do_mean_ctgs',1);
                        
            if plotmode == 6 || plotmode == 7
                warning('This still might not work properly');
                myfunc = @(x) (fcomplex2real(mean(x))); Nbs = 1000;
                group(i).data = bootstrp_and_scale(Nbs,myfunc,group(i).data')';
                group(i).datastats = bootstrp_and_scale(Nbs,myfunc,group(i).datastats')';
                clear myfunc Nbs
                
%                 % Bootstrap override
%                 group(i).data = fcomplex2real(mean(group(i).data,2));
%                 group(i).datastats = fcomplex2real(mean(group(i).datastats,2));
            
            end
            
%             if do_fisher
%                 group = group.myfunc (group, 'data', @fisher_trans);
%                 group = group.myfunc (group, 'datastats', @fisher_trans);
%             end

        case 3
            [group(i).data, group(i).funames] = get_grouped_data(group(i).ctgs,group(i).cells,ns,funames1D);
            timeband = get_stagesir(curr_stage_sfc)/1e3;
            [group(i).datastats, group(i).freqband_stats] = calc_pls_stats(t,group(i).data,timeband,'do_mean_ctgs',1);
            
        case 9
            [group(i).data, group(i).funames] = get_grouped_data(group(i).ctgs,group(i).cells,pls,funames1D);
            [group(i).datastats, group(i).freqband_stats] = calc_pls_stats(ftr,group(i).data,[freqband_stats],'do_mean_ctgs',1);
    end
    group(i).xdata=abscissa;
end

if ispaired
    group = grouppairs_merge(group,plotmode == 4 | plotmode == 7);
end

%% pls_stats
% % Calculate pls stat (note - we want this to go AFTER sorting section)
% % We also want to do this after loading data into groups, so we
% % don't duplicate any manipulation of pls (i.e. convert to complex->real).
pls_stats = calc_pls_stats(abscissa,pls(:,:,:),freqband_stats,'do_mean_ctgs',0);          % 11th column for all trials (need to make sure it's all good trials)
if plotmode == 7
    pls = fcomplex2real(pls);                % IF group is complex, convert it
    pls_stats = fcomplex2real(pls_stats);
end


if plot_SFC_spectra
    %% Plots Spectra

    figure;
    opts_PM3D = {'do_mean',domean,'do_zscore',do_zscore,'showErrorbars',showError};
    [h1, out.PM3Dcs] = plot_matrix3D_custstruct(abscissa,group,opts_PM3D,opts_PM3Dcs);
    
    
    switch plotmode
        %case 1; xlabel('Freq (Hz)'); xlim([0 160]); if do_zscore; ylabel('SFC - normalized'); else ylabel('SFC'); end; if sfc_mode >= 2.0 && sfc_mode <= 2.2; ylim([0.02 0.14]); end
        %case 2; xlabel('Freq (Hz)'); xlim([0 60]); ylabel('PSD');
        case 3; xlabel('Time (s)'); ylabel('Spikes/s'); add_stages_separators([],'k'); add_stages_names([],'k');
        %case {4,6}; xlabel('Freq (Hz)'); ylabel('Phase (rad)'); xlim([0 160]); 
        case 7; xlabel(['' num2str(freqband_corrcoef) 'Hz) coupling creq (Hz)']);
    end
end


if plot_SFC_spectra_stats
    %% Plot Stats
    fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
    h
    p
end

if fig_plot_histplot
    figure; [h1] = plot_histplot(group);
end

if fig_plot_venn
   figure; plot_venn_groups(group); 
end


if plot_spectrogram
    %% Plot Spectrogram
    for i = 1:length(group)
        [group(i).data, group(i).funames] = get_grouped_data(group(i).ctgs,group(i).cells,SPG,funames1D);
    end
    legend_arr = get_legendarr(group);
    
    % Load up data
            % Note: we need to do this outside the main plotting loop
            % in order to get gmin, gmax
    z = {group.data};
    for i = 1:length(z);
        if spectrogram_thresh > 0
            z{i} = z{i} > spectrogram_thresh;
            %z{i} = squeeze(mean(z{i},2)) * 100;   % Gives the % of cells above threshold; (Note: mean of a binary gives the percentage of 1's)
            z{i} = squeeze(sum(z{i},2)) * 1;
        else
            if floor(spect_mode) == 8
                z{i} = zscore(z{i},[],3);
            end 
            %z{i}=squeeze(mean(z{i},2));
            %z{i}=squeeze(median(z{i},2));
            z{i}=squeeze(quantile(z{i},0.75,2));      % Same as median
            
        end
    end
    [gmin, gmax] = get_plot_range(z);
    
    for i = 1:length(group)
        if isempty(group(i).data); continue; end
        figure;

        imagesc([min(T) max(T)],[min(F) max(F)], z{i});set(gca,'YDir','normal'); colorbar; %caxis([gmin,gmax]);
        add_stages_separators([],'w'); add_stages_names([],'w'); title(legend_arr{i})

        xlabel('Time (s)'); ylabel('Freq (Hz)');
    end
    
    % Plot animated for all cells
    if plot_spectrogram_individuals
        
        for i = 1:40; clf; imagesc([min(T) max(T)],[min(F) max(F)],squeeze(SPG(:,i,end,:))); set(gca,'YDir','normal'); add_stages_separators([],'w'); add_stages_names([],'w'); pause; end
        
        for i = 1:length(group)
            
            figure;        
            z=group(i).data;
            if spectrogram_thresh > 0
                z = z > spectrogram_thresh;
            end
            fprintf(['Group ' num2str(i) '\n']);
            
            if floor(spect_mode) == 3
                for j = 1:size(z,2);
                    clf;
                    imagesc([min(T) max(T)],[min(F) max(F)], squeeze(z(:,j,:)));set(gca,'YDir','normal'); colorbar;
                    add_stages_separators([],'w'); add_stages_names([],'w'); 
                    %title(['Group ' num2str(i) ' Cell #' num2str(j) '']);
                    pause
                end
            elseif floor(spect_mode) == 8
                for j = 1:size(z,2);
                    clf;
                    imagesc([min(T) max(T)],[min(F) max(F)], zscore(squeeze(z(:,j,:)),[],2));set(gca,'YDir','normal'); colorbar;
                    add_stages_separators([],'w'); add_stages_names([],'w'); 
                    %title(['Group ' num2str(i) ' Cell #' num2str(j) '']);
                    pause
                end
            end
            xlabel('Time (s)'); ylabel('Freq (Hz)');
            
            pause
        end
    end    
end

if plot_evoked
    %% Plot Evoked
    
    for i = 1:length(group)
        [group(i).data, group(i).funames] = get_grouped_data(group(i).ctgs,group(i).cells,ev,funames1D);
    end
    
    if plot_ev_raw
       
        figure;
        abscissa = 1:size(ev,1); abscissa = abscissa * get_dt - 2;
        opts_PM3D = {'do_mean',domean,'do_zscore',do_zscore,'showErrorbars',showError};
        [h1, out.PM3Dcs] = plot_matrix3D_custstruct(abscissa,group,opts_PM3D,opts_PM3Dcs);
        add_stages_separators([],'k'); add_stages_names([],'k'); 
        
    end
    
    if plot_ev_psd
        
        group_psd = group;
        for i = 1:length(group)
            group_psd(i).data = group(i).data((2600:3600),:);
            [group_psd(i).data, f] = psd_wrapper(group_psd(i).data,'fs',fs);
            group_psd(i).data = (group_psd(i).data);
        end
        opts_PM3D = {'do_mean',domean,'do_zscore',do_zscore,'showErrorbars',showError};
        [h1, out.PM3Dcs] = plot_matrix3D_custstruct(abscissa,group,opts_PM3D,opts_PM3Dcs);

    end
    
    if plot_ev_spect
        group_spect = group;
        for i = 1:length(group);
            [group_spect(i).data, F, T] = spect_wrapper(group(i).data,'fs',fs);
            group_spect(i).data = (group_spect(i).data);
            
            % figure; imagesc(T,F,zscore(group_spect(i).data,[],2)); set(gca,'YDir','normal'); colorbar
            figure; imagesc(T,F,log(group_spect(i).data)); set(gca,'YDir','normal'); colorbar
        end
        
    end
end

if plot_spectrogram_vs_ev
    

    ctgrange = 1:4;

    x = mean(ev(:,:,ctgrange),3); 
    y = mean(ns(:,:,ctgrange),3);
    z = squeeze(mean(SPG(:,:,ctgrange,:),3)); z = zscore(z,[],3);

    x = zscore(x);
    y = zscore(y);

    figure
    for i = 1:size(z,2)
        i
        clf;
        
        subplot(211); plott_fs(x(:,i),'b','fs',1/get_dt); hold on; plott_fs(y(:,i),'r','fs',1/get_dt); legend('lfp','spike rate')
        xlim([min(T) max(T)]); title(funames1D{i}); xlabel('Time(s)');
        subplot(212); imagesc([min(T) max(T)],[min(F) max(F)],squeeze(z(:,i,:))); set(gca,'YDir','normal');
        xlabel('Time(s)'); ylabel('Freq (Hz)');
        pause;
%         prompt = 'Type a comment or hit enter to continue q to quit: ';
%         returns{i} = input(prompt,'s');
%         if strcmp(returns{i},'q') || strcmp(returns{i},'Q'); break; end
        
    end


end

%% Plot Scatter
if scattergroups_sens_vs_sens
    %% %% Sens vs sens
    
    bad_scatter = bad_any | bad_zeros_sc;
    [group] = remove_bads_cells(bad_scatter,group);   % Remove additional data points that were associated with bad_zeros_sc (There should be none of these; see parse_sens.m)
    
    figure; plot(sens_sc(~bad_scatter,1),sens_sc(~bad_scatter,2),'k.');
%     figure; plott_fit(sens_sc(~bad_scatter,1),sens_sc(~bad_scatter,2),'k.');
    hold on; plot_scattergroups(sens_sc(:,1), sens_sc(:,2),group);
    hold on; plot([0 1.6],[0 1.6],'k')
    
    
%     xlim([-0.01 0.2]); ylim([-0.01 0.2])
    xlim([-0.00 0.9]); ylim([-0.00 0.9])
    xlabel(''); ylabel('')
    
    clear bad_scatter
end


if scattergroups_sens_vs_pls
    %% %% Pls vs sens
    mydata = pls_stats(:,end);
    bad_scatter = bad_any | bad_zeros_sc;
    [group] = remove_bads_cells(bad_scatter,group);   % Remove additional data points that were associated with bad_zeros_sc (There should be none of these; see parse_sens.m)
    

    if ~isreal(mydata)
        warning('Warning data not real. There is no point to use mode 6 for scatterplots, since there is no averaging - instead use mode 4. \n')
    end
    
    %figure; plot(sens_stat(~bad_scatter,sensitivity_coords1),mydata(~bad_scatter),'k.');
    figure; plott_fit(sens_sc(~bad_scatter,1),mydata(~bad_scatter),'k.');
    [R,P] = corrcoef(sens_sc(~bad_scatter,1),mydata(~bad_scatter));
    Rarr(i) = R(1,2); Parr(i) = P(1,2);
    hold on; plot_scattergroupdata(sens_sc(:,1) ,group, senscstr(1));
    
end

%% Plot Regress
if regress_pls_vs_sens
    bad_regress = bad_any | bad_zeros_reg;

    if ~isreal(sens_reg)
        warning('Warning data not real. There is no point to use mode 6 for scatterplots, since there is no averaging - instead use mode 4. \n')
    end

    reg_struct = plots_regress(sens_reg(~bad_regress,:),sensregnames(:),group,'.');
    %sens_reg = [sens_reg(:,1:end-1) score(:,4) sens_reg(:,end)];
    %sensregstr = {sensregstr{1:end-1} 'Score' sensregstr{end}};
    %plots_regress(sens_reg(~bad_regress,:),sensregstr(:),group,'.');
    temp=repmat({freqband_stats},1,length(reg_struct));
    [reg_struct.freqband_stats] = temp{:};
end

if plot_fisher_test_suite
     bad_regress = bad_any | bad_zeros_reg;
     if sum(bad_zeros_reg) > 0; warning('Bad zeros regress >0'); keyboard; end
     close all
     X0 = sens_reg(~bad_regress,end);
     X0=X0(:);
     test_fisher_suite
end

%% Exploratory Plots
if find_archetypal_beta_cells
 
    plot_archetypal_betas(group,pls,spc);   % Note - must be in group = gr_all mode, so length(group) = 16
    % Do stuff manually
    n = find(any(group(9).cells,2));
    [fnum, unum, curr_funame] = autorun_unum2fnum(n(1))
    [fnum, unum, curr_funame] = autorun_unum2fnum(n(2))
    [fnum, unum, curr_funame] = autorun_unum2fnum(n(3))
    figure; plott_ani(squeeze(pls(:,n,:)))
    
end




if plot_find_high_beta_cells
    for i = 1:length(group)
        
        ds = get_freq_stats(group(i).data,f,freqband_stats);
        figure; hist(ds,50);
        xlabel('Beta SFC'); ylabel('# cells');
        
        
        indhigh = find_high_beta_cells(group(i),f,0.08,freqband_stats);
        figure; plot(group(i).data(:,indhigh)); xlim([0 160]);
        xlabel('F (Hz)'); ylabel('SFC'); legend('High SFC cells');
        
    end
end

if phi_survey_ctgs
    %% Phi Survey
    
    [ph1a, ph2a, delphi] = get_phasestats_ctgs(f,Cave,phiave,freqband_stats, [9 10],group);
    
    chosen_cells = ~bad_any;
    ph_senscoords=10;
    figure; plott_fit(sens(chosen_cells,ph_senscoords),abs(delphi(chosen_cells)),'k.');
    hold on; plot_scattergroups(sens(:,ph_senscoords), abs(delphi),group);
    
%     figure; plot(Ca1(chosen_cells),Ca2(chosen_cells),'k.');
%     hold on; plot_scattergroups(Ca1,Ca2,group(1:2));
    
    figure; plot(ph1a(chosen_cells),ph2a(chosen_cells),'k.');
    hold on; plot_scattergroups(ph1a,ph2a,group(1:2));
    
%     figure; plot3(ph1a(chosen_cells),ph2a(chosen_cells),Ca(chosen_cells),'k.');
   
    %% Regression
    
%     ph_regress_coords=[1,2,3,4,8,12];
%     figure; chosen_cells = chosen_cells(:);
%     plots_regress(sens(chosen_cells,ph_regress_coords),ph3a(chosen_cells),sensnames(ph_regress_coords),group,'.',src(chosen_cells,end));
    
    
end


if survey_boundary_vs_categorization_index
    %% Phi Survey
    
    if ~isempty(strfind(fname_suffix_sfc,'bndry'))

        % Chosen cells
        chosen_cells = true(1,size(bad_any,1)) & ~bad_any;     % All
    %     chosen_cells = ismonkeyL & ~bad_any;               % Monkey L
    %     chosen_cells = ~ismonkeyL & ~bad_any;              % Monkey O

        % Load phase differences for SFC
        ctgs = [5,6,5];
        [ph1a, ph2a, delphi] = get_phasestats_ctgs(f,Cave,phiave,freqband_stats, ctgs,group);

        % Phi plot
        figure; plot(ph1a(chosen_cells),ph2a(chosen_cells),'k.');
        hold on; plot_scattergroups(ph1a,ph2a,group(1:2));
        xlabel(['Phi Ctg ' num2str(ctgs(1))]); ylabel(['Phi Ctg ' num2str(ctgs(2))]);
        title('Phase in different schemes');


        % Load categorization index
        %Sch A Rel
        ind = 9:9+5;
        catFR = statsFR(:,ind); catFR_names = statsnames(ind);

        WCD = [catFR(:,1)-catFR(:,2) catFR(:,2)-catFR(:,3) catFR(:,4)-catFR(:,5) catFR(:,5)-catFR(:,6)];
        WCD = abs(WCD);
        WCD = mean(WCD,2);

        BCD = mean(abs(catFR(:,3)-catFR(:,4)),2);

        CI_SchARel = (BCD-WCD) ./ (BCD+WCD);


        % Sch B Rel
        ind = 21:21+5;
        catFR = statsFR(:,ind); catFR_names = statsnames(ind);

        WCD = [catFR(:,1)-catFR(:,2) catFR(:,2)-catFR(:,3) catFR(:,4)-catFR(:,5) catFR(:,5)-catFR(:,6)];
        WCD = abs(WCD);
        WCD = mean(WCD,2);

        BCD = mean(abs(catFR(:,3)-catFR(:,4)),2);

        CI_SchBRel = (BCD-WCD) ./ (BCD+WCD);

        CI_mean = (CI_SchARel+CI_SchARel)/2;

        figure; plott_fit(CI_mean(chosen_cells),delphi(chosen_cells),'.'); xlabel('Category Index'); ylabel('Phase Difference');
    % % %     % delphi2=  ph1a-ph2a;
    % % %     % figure; plott_fit(CI_mean(chosen_cells),delphi2(chosen_cells),'.');


        [Ca] = calc_pls_stats(f,pls,freqband_stats,'do_mean_ctgs',0);
        Cadiff = abs(Ca(:,ctgs(2)) - Ca(:,ctgs(1)));
        figure; plott_fit(CI_mean(chosen_cells),Cadiff(chosen_cells),'.');  xlabel('Category Index'); ylabel('Coherence Difference');
    % % %     figure; plott_fit(CI_mean(chosen_cells),Ca(chosen_cells,end),'.');  xlabel('Category Index'); ylabel('Coherence');
    else
        warning('Should be in bndry mode for this to work (sfc_mode=2.2415101;)');
    end
end


if Cave_survey_ctgs
%% Cave_survey    
    [Ca] = calc_pls_stats(f,pls,freqband_stats,'do_mean_ctgs',0);
    Ca1 = mean(Ca(:,[9]),2);
    Ca2 = mean(Ca(:,[10]),2);
    chosen_cells = ~bad_any;
    
    figure; plott_fit(Ca1(chosen_cells),Ca2(chosen_cells),'k.');
    
    
    [R,P] = corrcoef(Ca1(chosen_cells),Ca2(chosen_cells));
    Rarr(i) = R(1,2); Parr(i) = P(1,2);

    hold on; plot_scattergroups(Ca1,Ca2,group(1:2));
    xlabel('Beta SFC SchA');ylabel('Beta SFC SchB');
    
    [R,P] = corrcoef(Ca);
    P = P - eye(size(P));
    %figure; imagesc(R); set(gca,'YDir','normal'); colorbar; caxis([0 1]); xlabel('Ctgs'); ylabel('Ctgs'); title(['Corr Coef with p < ' num2str(max(P(:)))]);
    %figure; imagesc(P); set(gca,'YDir','normal'); colorbar; caxis([0 1]); xlabel('Ctgs'); ylabel('Ctgs');
    
    
end


if plot_Sch_sens
    
    if curr_stage_sfc == 3
        s_boot = wrkspc_buffer.(['boot_' mode2modename(boot_mode)]).(stagename(curr_stage_sfc));
    else
        s_boot = [];    % Since we don't have any stage2 data
    end
    [out.PPS] = plot_permute_scatter(...
        wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
        s_boot, ...
        bad_any,freqband_stats, opts_PPS);
end


if plot_PCA
    %% PCA
    %%%%%%% Get phases
  
    [ph1a, ph2a, delphi, Ca] = get_phasestats_ctgs(f,Cave,phiave,freqband_stats, [9 10],group);
    mydata = pls_stats(:,end);
    ph = delphi;
    
    %delphi = abs(delphi);
    figl; subplotcols(2,1); hist(delphi,20); xlabel('Phase (rad)'); ylabel('NCells')
    subplotcols(2,2); plot(delphi(~bad_any),Ca(~bad_any),'k.')
    %hold on; plot_scattergroupdata(delphi ,group, 'Phase (rad)');
    
    hold on; plot_scattergroups(delphi, Ca,group(1:2));
    
    %%%%%%%%

    %X = [pls_stats(:),sens, sens_alt];
    %X = [pls_stats(:) sens(:,[1:4,9]), sens_alt(:,[1:4,9])];
 
%     PCA_coords = [1:4,9,12];
%     X = [(sens(:,PCA_coords)) (sens_alt(:,PCA_coords))]; Xnames = {sensnames{PCA_coords}}; Xnames_alt = {sensnames_alt{PCA_coords}};
%     X = [mydata, delphi, X]; Xnames = {'Pls','Phi', Xnames{:}};           % Add phi values
%     Xnames_full = {Xnames{:} Xnames_alt{:}}
    
    reorder = 1;
    PCA_coords = [1:4,10];
%     PCA_coords = [1:2, 10];
    pca_stages = [-2 -1 0 2 3];
    clear srpca sspca senspca sensnamespca
    for i = 1:length(pca_stages)
        [~,srpca{i},sspca{i},pr_fnames,pr_funames1D] = extract_spsr(pca_stages(i),bad_any);
        [senspca{i}, sensnamespca{i}] = calc_sens(srpca{i},sspca{i},pca_stages(i),use_morph_sens,sens_std_normalize);
        sensnamespca{i} = strcat(['Stg' num2str(pca_stages(i)) ':'], sensnamespca{i}); 
    end
    Xnames = vertcat(sensnamespca{:}); Xnames = Xnames(:,PCA_coords); Xnames = permute(Xnames,[3 2 1]);
    X_gr = cat(3,senspca{:}); X_gr = X_gr(:,PCA_coords,:);
    if reorder; X_gr = permute(X_gr,[1 3 2]); Xnames = permute(Xnames,[1 3 2]); end
    sz=size(X_gr);
    X_1D = reshape(X_gr, [sz(1),sz(2)*sz(3)]);
    Xnames_1D = reshape(Xnames,[1,sz(2)*sz(3)]);
    
    figure;
    N = size(X_gr,3);
    
    
    % Get unique portion of "rows" from Xnames (2nd dimension)
    temp = Xnames(1,1,:);
    temp = cellfunu(@(s) strsplit(s,':'), temp); temp=vertcat(temp{:});  %Split names up into stage and sensnames parts
    uniquecol = [find(sum(strcmp(temp{1,1},temp)) == 1) find(sum(strcmp(temp{1,2},temp)) == 1)];
    Xnames_rows = temp(:,uniquecol);
    
    
    %hsp=subplot_grid(N,1);
    for i = 1:N
        %hsp.set_gca(i);
        subplot(N+1,1,i);bar_matrix3D(squeeze(X_gr(:,:,i)),'active_dim',1);
        legend(Xnames_rows{i});
        set(gca,'XTickLabel',[])
    end
    subplot(N+1,1,N+1); bar(zeros(5,1))
    
    % Get unique "columns" from Xname (3rd dimension)
    temp = Xnames(1,:,1);
    temp = cellfunu(@(s) strsplit(s,':'), temp); temp=vertcat(temp{:});  %Split names up into stage and sensnames parts
    uniquecol = [find(sum(strcmp(temp{1,1},temp)) == 1) find(sum(strcmp(temp{1,2},temp)) == 1)];
    Xnames_cols = temp(:,uniquecol);
    
    %Xnames_cols = strrep(Xnames_cols,'SchAR','');
    xticklabel_rotate(1:size(Xnames,2),90,Xnames_cols)
    clear temp uniquecol Xnames_cols
    
%     for i = 1:
    
%     X=[]; Xnames_full={};
%     for i = 1:length(pca_stages)
%         X = [X senspca{i}(:,PCA_coords)];
%         Xnames_full = {Xnames_full{:} sensnamespca{i}{PCA_coords}};
%     end
    
    
%     index=1:549;
    index = ~bad_any(:);
%     index = ~bad_any(:) & (any([group(1).cells,group(2).cells],2));
%     index = ~bad_any(:) & Ca > 0.00;
    X_1D = X_1D(index,:);
    %X = [score_ph(:,3), X];
    %X = log(X);
    zX = zscore((X_1D)); 
    
    [coeff,score,latent] = princomp(zX);
%     
    figl; set(gcf,'Position',[1 26 1600 869]);
    hsp=subplot_grid(4,1);
    hsp.set_gca(1); bar(coeff(:,1)); set(gca,'XTickLabel',[])
    
    hsp.set_gca(2); bar(coeff(:,2)); set(gca,'XTickLabel',[])
    hsp.set_gca(3); bar(coeff(:,3)); set(gca,'XTickLabel',[])
    %hsp.set_gca(4); bar(zeros(size(coeff(:,1)))); xticklabel_rotate(1:length(Xnames_1D),90,Xnames_1D);


%     figl; set(gcf,'Position',[ 1          26        1600         869]);
%     subplotrows(3,1); bar(coeff(:,1)); xticklabel_rotate(1:length(Xnames_full),90,Xnames_full)
%     subplotrows(3,2); bar(coeff(:,2)); subplotrows(3,3); bar(coeff(:,3));
   

    
    figure; plot(score(:,1),score(:,2),'k.');
    hold on; plot_scattergroups(score(:,1), score(:,2),group,index);
    %hold on; plot([0 20],[0 20],'k');
    
    figure; plot3(score(:,1),score(:,2),score(:,3),'k.');
    
    
    figure; bar(mean(X_1D));
    xticklabel_rotate(1:length(Xnames_1D),90,Xnames_1D)
    %set(gca,'YTick',1:length(Xnames_full));
    %set(gca,'YTickLabel',Xnames_full);
    
    [R] = corrcoef(X_1D);
    figure; imagesc(R);
    xticklabel_rotate(1:length(Xnames_1D),90,Xnames_1D)
    set(gca,'YTick',1:length(Xnames_1D));
    set(gca,'YTickLabel',Xnames_1D);
    
    %%
    %figure; chosen_cells = ~bad_any(:);
    %plots_regress(X(:,2:end),X(:,1),{Xnames{2:end},Xnames_alt{1:end}},group,'.',src(chosen_cells,end));
    
    
end
%%

%% A bunch of random tests
if 0
    keyboard
    %% Analyze specialization indices
    statsnames = stats.stage2.stats(1).m.ctgsetnames;
    statsFR = stats.stage2.mu_arr;
    statsFR = squeeze(statsFR);
    
    chosen_cells = ones(size(ismonkeyL));     % All
%     chosen_cells = ismonkeyL;               % Monkey L
%     chosen_cells = ~ismonkeyL;              % Monkey O
    chosen_cells = chosen_cells(:);
    

    
    % Specialization Indices - A number between 0 and 1 that tells how
        % specialzed a cell is
    clc
    Sx = abs(statsFR(:,1)-statsFR(:,3)).^1; Sy = abs(statsFR(:,2)-statsFR(:,4)).^1; speciARel = abs(Sx - Sy) ./ (Sx + Sy);
    Sx = abs(statsFR(:,5)-statsFR(:,6)).^1; Sy = abs(statsFR(:,7)-statsFR(:,8)).^1; speciBRel = abs(Sx - Sy) ./ (Sx + Sy);
    
    Sx = abs(statsFR(:,1)-statsFR(:,2)).^1; Sy = abs(statsFR(:,3)-statsFR(:,4)).^1; speciANR = abs(Sx - Sy) ./ (Sx + Sy);
    Sx = abs(statsFR(:,5)-statsFR(:,7)).^1; Sy = abs(statsFR(:,6)-statsFR(:,8)).^1; speciBNR = abs(Sx - Sy) ./ (Sx + Sy);
    SpecI=[speciARel speciBRel  speciANR speciBNR];         % Pack into one array
    
    % Specialization Rankings - Ranks morph quandrants in order of firing rate
        % Sch A
    X = [statsFR(:,1) statsFR(:,2) statsFR(:,3),statsFR(:,4)];
    X = sort(X,2,'descend');
    X = X ./ repmat(X(:,1),1,size(X,2));
        % Sch B
    Y = [statsFR(:,5) statsFR(:,7) statsFR(:,6),statsFR(:,8)];
    Y = sort(Y,2,'descend');
    Y = Y ./ repmat(Y(:,1),1,size(Y,2));
        % Sch A & B
    Z = [statsFR(:,1) statsFR(:,2) statsFR(:,3),statsFR(:,4),statsFR(:,5) statsFR(:,7) statsFR(:,6),statsFR(:,8)];
    Z = sort(Z,2,'descend');
    Z = Z ./ repmat(Z(:,1),1,size(Z,2));
    
    %% Specialization Indices vs Groups
    indA = (speciARel(group(3).cells & chosen_cells));
    indB = (speciBRel(group(2).cells & chosen_cells));
    indC = (speciANR(group(4).cells & chosen_cells));
    indD = (speciBNR(group(4).cells & chosen_cells));
    
    quant_thresh = 0.7;
    indA = speciARel(sens(chosen_cells,1) > quantile(sens(chosen_cells,1),quant_thresh));
    indB = speciBRel(sens(chosen_cells,2) > quantile(sens(chosen_cells,2),quant_thresh));
    indC = speciANR(sens(chosen_cells,1) <= quantile(sens(chosen_cells,1),quant_thresh));
    indD = speciBNR(sens(chosen_cells,2) <= quantile(sens(chosen_cells,2),quant_thresh));

    figure; hist(indA)
    figure; hist(indB)
    
    mA = median(indA)
    mB = median(indB)
    mC = median(indC)
    mD = median(indD)

    % Stats
    [p] = ranksum(indA,indB)
    [p] = ranksum(indB,indD)
    
    % Sample plots
    if 0
        i=548; Xsamp = [statsFR(i,1) statsFR(i,3); statsFR(i,2) statsFR(i,4)];
        figure; imagesc(Xsamp)
        speciARel(i)

        i=43; X = [statsFR(i,5) statsFR(i,7); statsFR(i,6) statsFR(i,8)];
        figure; imagesc(X)
        speciBRel(i)
    end
    
    % PCA on specialization index
    if 0
        % Not much interesting here; just 
        [coeff score latent] = princomp(zscore(SpecI));
    end
    
    %% PCA on Firing rates
    [coeff score latent] = princomp(zscore(statsFR(:,1:8)));
    Xregress = [abs(score(:,1:5)) pls_stats'];
    figure; bar(coeff(:,4)); % Diagonally sensitive neurons!
    plots_regress(Xregress(~bad_any,:),{'FR','A1vA2','B1vB2','Diag','SchAB','Pls'});
    
    

    %% Specialization Rankings Clustering
    
    %figure; plot( (mean(X(:,1:2),2) - mean(X(:,3),2)), (mean(X(:,1:3),2) - mean(X(:,4),2)),'.')
    
    figure; plot(X(:,1)-X(:,2),X(:,2)-X(:,3),'.')
    %figure; plot(X(:,2),X(:,3),'.')

    
    %% Specialization Rankings vs Groupings
    clc
    
    % Cell grouping method 1 - Groups
    indA = 1-X(group(3).cells & chosen_cells,2);
    indB = 1-Y(group(2).cells & chosen_cells,2);
    indC = 1-X(group(4).cells & chosen_cells,2);
    indD = 1-Y(group(4).cells & chosen_cells,2);
    
    % Cell grouping method 2 - Raw Thresholds
%     thresh = 0.15;
%     indA = 1-X(sens(chosen_cells,1) > thresh ,2);
%     indB = 1-Y(sens(chosen_cells,2) > thresh ,2);
    
    % Cell grouping method 3 - Quantiles
    quant_thresh = 0.9;
    indA = 1-X(sens(chosen_cells,1) > quantile(sens(chosen_cells,1),quant_thresh) ,2);
    indB = 1-Y(sens(chosen_cells,2) > quantile(sens(chosen_cells,2),quant_thresh) ,2);
    indC = 1-X(sens(chosen_cells,1) <= quantile(sens(chosen_cells,1),quant_thresh) ,2);
    indD = 1-Y(sens(chosen_cells,2) <= quantile(sens(chosen_cells,2),quant_thresh) ,2);
    
    figure; hist(indA)
    figure; hist(indB)
%     
    mA = median(indA)
    mB = median(indB)
    mC = median(indC)
    mD = median(indD)

    [p h] = ranksum(indA,indB)
    
    %% Specialization index vs monkey
    
        % Sch A
    quant_thresh = 0.75;
    
    % All cells
    indA = speciARel(ismonkeyL);
    indB = speciARel(~ismonkeyL);
    
    % Specialist cells
    indA = speciARel(ismonkeyL' & sens(:,1) > quantile(sens(ismonkeyL,1),quant_thresh));
    indB = speciARel(~ismonkeyL' & sens(:,1) > quantile(sens(~ismonkeyL,1),quant_thresh));
    
    
    figure; hist(indA)
    figure; hist(indB)
%     
    mA = median(indA)
    mB = median(indB)

    [p h] = ranksum(indA,indB)
    
        % Sch B
    prompt = 'Hit enter to continue or q to quit: '; ret = input(prompt,'s'); if strcmp(ret,'q') || strcmp(ret,'Q'); error('Quit by user'); end
    
    % All cells
    indA = speciBRel(ismonkeyL);
    indB = speciBRel(~ismonkeyL);
    
    % Specialist cells
    indA = speciBRel(ismonkeyL' & sens(:,2) > quantile(sens(ismonkeyL,2),quant_thresh));
    indB = speciBRel(~ismonkeyL' & sens(:,2) > quantile(sens(~ismonkeyL,2),quant_thresh));
    
    figure; hist(indA)
    figure; hist(indB)
%     
    mA = median(indA)
    mB = median(indB)

    [p h] = ranksum(indA,indB)
    
    
    %% Specialization ranking vs monkey
    
        % Sch A
    quant_thresh = 0.5;
    indA = 1-X(ismonkeyL,2);
    indB = 1-X(~ismonkeyL,2);
    
%     indA = 1-X(ismonkeyL' & sens(:,2) > quantile(sens(ismonkeyL,2),quant_thresh) & sens(:,1) > quantile(sens(ismonkeyL,1),quant_thresh),2);
%     indB = 1-X(~ismonkeyL' & sens(:,2) > quantile(sens(~ismonkeyL,2),quant_thresh) & sens(:,1) > quantile(sens(~ismonkeyL,1),quant_thresh),2);
    
    
    figure; hist(indA)
    figure; hist(indB)
%     
    mA = median(indA)
    mB = median(indB)

    [p h] = ranksum(indA,indB)
    
        % Sch B
    prompt = 'Hit enter to continue or q to quit: '; ret = input(prompt,'s'); if strcmp(ret,'q') || strcmp(ret,'Q'); error('Quit by user'); end
    
    indA = 1-Y(ismonkeyL,2);
    indB = 1-Y(~ismonkeyL,2);
    
    figure; hist(indA)
    figure; hist(indB)
%     
    mA = median(indA)
    mB = median(indB)

    [p h] = ranksum(indA,indB)
    
    %% Specialization Rankings Regression
    %ch = true(size(sens,1),1);
    ch = group(2).cells | group(3).cells | group(1).cells;
    plots_regress([sens(ch,1),sens(ch,2),X(ch,2)],{'SensA','SensB','Speicalization'});
    %figure; plott_fit(sens(ch,1),X(ch,2),'.');
    

    
%     could these metrics still be biased by changes in firing rate??


    %% Specialization Rankings vs PLS - For sample neurons, neurons specializing in 2 categories neurons are sig; for delay "2 and 3 specialists are sig
    ch = true(size(sens,1),1);
%     ch = group(2).cells | group(3).cells | group(1).cells;
    X_temp = X;
    plots_regress([X_temp(ch,1)-X_temp(ch,2) X_temp(ch,2)-X_temp(ch,3) X_temp(ch,3)-X_temp(ch,4) sens_alt(ch,9) pls_stats(ch)' ],{'Spec1','Spec2','Spec3','FR','PLS'});
    
    %figure;plot(X_temp(ch,2),X_temp(ch,3),'.')
    %figure; plot(X_temp(ch,1)-X_temp(ch,2), X_temp(ch,2)-X_temp(ch,3),'.')
    
end


if 0
    keyboard
    %% Categorization index for A vs B
    
    chosen_cells = true(size(ismonkeyL));     % All
%     chosen_cells = ismonkeyL;                 % Monkey L
%     chosen_cells = ~ismonkeyL;              % Monkey O
    chosen_cells = chosen_cells(:);
    
    % Categorization index
    chosen_cells = true(1,size(choup,1)) & ~bad_any;     % All
%     chosen_cells = ismonkeyL & ~bad_any;               % Monkey L
%     chosen_cells = ~ismonkeyL & ~bad_any;              % Monkey O
    
    %Sch A Rel
    ind = 9:9+5;
    catFR = statsFR(:,ind); catFR_names = statsnames(ind);
    
    WCD = [catFR(:,1)-catFR(:,2) catFR(:,2)-catFR(:,3) catFR(:,4)-catFR(:,5) catFR(:,5)-catFR(:,6)];
    WCD = abs(WCD);
    WCD = mean(WCD,2);
    
    BCD = mean(abs(catFR(:,3)-catFR(:,4)),2);
    
    CI_SchARel = (BCD-WCD) ./ (BCD+WCD);
    
    
    % Sch B Rel
    ind = 21:21+5;
    catFR = statsFR(:,ind); catFR_names = statsnames(ind);
    
    WCD = [catFR(:,1)-catFR(:,2) catFR(:,2)-catFR(:,3) catFR(:,4)-catFR(:,5) catFR(:,5)-catFR(:,6)];
    WCD = abs(WCD);
    WCD = mean(WCD,2);
    
    BCD = mean(abs(catFR(:,3)-catFR(:,4)),2);
    
    CI_SchBRel = (BCD-WCD) ./ (BCD+WCD);
    
    
    indA = (CI_SchARel(group(3).cells & chosen_cells));
    indB = (CI_SchBRel(group(2).cells & chosen_cells));
    
    quant_thresh = 0.7;
    indA = CI_SchARel(sens(chosen_cells,1) > quantile(sens(chosen_cells,1),quant_thresh));
    %indB = CI_SchARel(sens(chosen_cells,1) < quantile(sens(chosen_cells,1),0.5));
    indB = CI_SchBRel(sens(chosen_cells,2) > quantile(sens(chosen_cells,2),quant_thresh));
    
    
    figure; hist(indA)
    figure; hist(indB)
    
    mA = median(indA)
    mB = median(indB)

    % Stats
    [p h] = ranksum(indA,indB)
    
    
    %% Categorization index versus monkeys
    quant_thresh=0.75;
    
        % Sch A
        
    indA = (CI_SchARel(ismonkeyL' & group(3).cells));
    indB = (CI_SchARel(~ismonkeyL' & group(3).cells));
    
%     indA = (CI_SchARel(ismonkeyL' & sens(:,1) > quantile(sens(ismonkeyL,1),quant_thresh)));
%     indB = (CI_SchARel(~ismonkeyL'  & sens(:,1) > quantile(sens(~ismonkeyL,1),quant_thresh)));
    
    figure; hist(indA)
    figure; hist(indB)
    
    mA = median(indA)
    mB = median(indB)

    % Stats
    [p h] = ranksum(indA,indB)
    prompt = 'Hit enter to continue or q to quit: '; ret = input(prompt,'s'); if strcmp(ret,'q') || strcmp(ret,'Q'); error('Quit by user'); end
    
        % Sch B
    indA = (CI_SchBRel(ismonkeyL' & group(2).cells));
    indB = (CI_SchBRel(~ismonkeyL' & group(2).cells));    
%     indA = (CI_SchBRel(ismonkeyL' & sens(:,2) > quantile(sens(ismonkeyL,2),quant_thresh)));
%     indB = (CI_SchBRel(~ismonkeyL' & sens(:,2) > quantile(sens(~ismonkeyL,2),quant_thresh)));
    
    figure; hist(indA)
    figure; hist(indB)
    
    mA = median(indA)
    mB = median(indB)

    % Stats
    [p h] = ranksum(indA,indB)
    
end


if 0
    %% Input vs Memory neurons -- Sample and Memory neurons have more overlap for A than for B
    
    indA = sens_alt(group(3).cells,1);      % Delay sensitivity to A
    indB = sens_alt(group(2).cells,2);      % Delay sensitivity to B
    %
%     thresh=0.15;
%     indA = sens_alt(sens(:,1) > thresh ,1);
%     indB = sens_alt(sens(:,2) > thresh ,2);
    
    quant_thresh=0.9;
    indA = sens_alt(sens(:,1) > quantile(sens(:,1),quant_thresh) ,1);
    indB = sens_alt(sens(:,2) > quantile(sens(:,2),quant_thresh) ,2);

    
    figure; hist(indA)
    figure; hist(indB)
    
    mA = median(indA)
    mB = median(indB)

    % Stats
    [p h] = ranksum(indA,indB)
    
    
end

if 0
    %% Compare how "useful" A and B neurons are in their non-preferred schemes
        % No real difference!
    
    indA = sens(group(3).cells,2);  % Scheme B sensitivity of Scheme A cells
    indB = sens(group(2).cells,1);  % Scheme A sensitivity of Scheme B cells
    indC = sens(group(4).cells,1);
    indD = sens(group(4).cells,2);
    
    figure; hist(indA)
    figure; hist(indB)
    
    mA = median(indA)
    mB = median(indB)
    
    mC = median(indC)
    mD = median(indD)

    % Stats
    [p h] = ranksum(indA,indB)
    
    
end
    

if 0
    keyboard
    %% Test changes in alpha and beta sample and delay
    CaveS = wrkspc_buffer.sfc_2_20151.stage2.Cave(:,:,end);
    CaveD = wrkspc_buffer.sfc_2_20151.stage3.Cave(:,:,end);
    Salpha = calc_pls_stats(f,CaveS,[10 12],'do_mean_ctgs',0,'do_mean_freq',0);
    Sbeta = calc_pls_stats(f,CaveS,freqband_stats,'do_mean_ctgs',0,'do_mean_freq',0);
    
    Dalpha = calc_pls_stats(f,CaveD,[10 12],'do_mean_ctgs',0,'do_mean_freq',0);
    Dbeta = calc_pls_stats(f,CaveD,freqband_stats,'do_mean_ctgs',0,'do_mean_freq',0);
    
    ismonkeyL_temp=~cellfun(@isempty,(strfind(wrkspc_buffer.sfc_2_20151.stage3.funames1D,'L')));
    
    %% Monkey L
    ind = ismonkeyL_temp;
    SAc = Salpha(ind);
    SBc = Sbeta(ind);
    DAc = Dalpha(ind);
    DBc = Dbeta(ind);
    [p h] = signrank(SAc,DAc)
    [p h] = signrank(SBc,DBc)
    figure;
    subplot(121); bar_matrix3D([SAc', DAc']); title('Alpha')
    subplot(122); bar_matrix3D([SBc', DBc']); title('Beta');
    
    %% Monkey O
    ind = ~ismonkeyL_temp;
    SAc = Salpha(ind);
    SBc = Sbeta(ind);
    DAc = Dalpha(ind);
    DBc = Dbeta(ind);
    [p h] = signrank(SAc,DAc)
    [p h] = signrank(SBc,DBc)
    figure;
    subplot(121); bar_matrix3D([SAc', DAc']); title('Alpha')
    subplot(122); bar_matrix3D([SBc', DBc']); title('Beta');
    
    %% Monkey all
    ind = true(size(ismonkeyL_temp));
    SAc = Salpha(ind);
    SBc = Sbeta(ind);
    DAc = Dalpha(ind);
    DBc = Dbeta(ind);
    [p h] = signrank(SAc,DAc)
    [p h] = signrank(SBc,DBc)
    figure;
    subplot(121); bar_matrix3D([SAc', DAc']); title('Alpha')
    subplot(122); bar_matrix3D([SBc', DBc']); title('Beta');
    
end


if bootstrap_validation_KStest
    %% Analyze bootstrap tests
    vars_pull(wrkspc_buffer.(buff_fieldname).(stagename(3)));
    ksp = horzcat(ks.p);
    ksstat = horzcat(ks.ksstat);
    ztp = horzcat(zerotest.p);
    ztpd = horzcat(zerotest.pd);
    sig = dCave_bs_sig;
    
    %% Plot KS Stat
    ind = find(f >= 30,1,'first');
    
    figure; subplot(121); imagescy(ksstat);  title('Ksstat'); ylabel('Freq'); xlabel('Cell #'); colorbar; subplot(122);imagescy(ksp < 0.05); ylabel('Freq'); xlabel('Cell #');title('KS Test (non-Gaussian)');colorbar;
    figure; subplot(121); plot(ksstat(ind,:)); title('Ksstat @ 24-36 Hz '); ylabel('ksstat'); xlabel('Cell #'); subplot(122); plot(sum(ksp <0.05,2)); title('KS Test'); xlabel('Freq(Hz)'); yabel('# cells H1'); 
    
    %% Plot P Difference
    pd = 100;    % percent
    x = abs(Cave1 - Cave2);
    
    % Implement method to get p value
    
    mu = (Cave1 + Cave2) / 2;
    delta = mu * (pd/100);                % Reject null hypothesis for changes in Cave less than this value
    
    sz = size(x);
    p = 1 - normcdf(delta - x,zeros(sz),sig);
    h = p < 0.05;
    
    figure; subplot(121); imagescy(p);  title(['p > ' num2str(pd) '% change']); ylabel('Freq'); xlabel('Cell #'); colorbar; subplot(122);imagescy(h); ylabel('Freq'); xlabel('Cell #');title(['h, alpha=0.05']); colorbar;
    figure; plot(sum(h,2)); ylabel('# cells no change'); xlabel('Freq(Hz)'); title(['h, alpha=0.05']);  
    %%
    keyboard 
end

% Package data if this script is called by a parent function
if is_calledby
    out.wrkspc_buffer = wrkspc_buffer;
    out.group = group;
    if exist('reg_struct','var'); out.reg_struct=reg_struct; end
    out.bad_any = bad_any;
    out.pls = pls;
    out.pls_stats = pls_stats;
end


% Clear everything except wrkspc
clearvars -except wrkspc_buffer out group gr f ftr pls fv % Cannot comment this out, because variables will not refresh

