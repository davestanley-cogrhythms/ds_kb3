% This code works with the results of permutation tests
% Tests differences in Amp & Phase across category schemes.
% Data should be generated by run_analysis.m AND with:
%     sfc_mode = *.?????1
%     OR, equivalently, folders ending in _permu

%% Stage
% clear
clc
run_setdefaultfig
addpath(genpath('./funcs_plots_preferred/'));

% Data loading parameters
recalc_md = 0;
recalc_sfc = 1;
recalc_permute = 0;
recalc_bootstrap = 0;
reload_data = 0;
if get_iscromer; stage_range = [2,3];
else stage_range = [2,3];
end
sfc_mode=22.4918103;
perm_mode=22.4614111;
boot_mode = 22.401312;
    isolate_clustered_protocols = 1; % 0 for all, 1 for cells with clustered switch trials, 2 for only cells with no clustering of switch trials
spect_mode = 3.20101;

%perm_mode=sfc_mode+0.000001;   % Make sure these are aligned

gr_submode = [1];             % Which of the ctgs to plot from the bndry, RT modes (i.e. congruent/incongruent; 50 or 60% morphs

[~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
fprintf(['Chosen SFC mode: ' fname_suffix_sfc ' \n']);

% Load preferred units
    % Stage selection
    curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;

    % Unit exclusion
    exclude_clipping = 1;
    exclude_60 = 0;
    exclude_nans = 1;
    excludeL = 0; 
    excludeO = 0; 
    
    % Pls switches
    do_fisher=0;
%     freqband_stats = [1 4];
%     freqband_stats = [10 12];
    freqband_stats = [16 20];
    plotmode = 1;               % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    perm2pls = 1;               % Instead of showing raw SFC values, show % cells successfully passing permutation test
        perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
        perm2pls_dophi = 0;
    collapse_pls_to_days = 0;
    
    opts_PM3Dcs.paperfig_mode=0;
    opts_PM3Dcs.do_sgolay=0;
    opts_PSC.paperfig_mode=0;
    opts_PPS.paperfig_mode = 0;
    opts_PPS.sens_chosen_condition=1;
    opts_PPS.hide_boot_results=1;
    opts_PPSC.paperfig_mode=0;

%% Plotting parameters
    % General settings


% Sorting of data (this is kind of slow)
sort_on = 0;


%% Plot Switches

% General
domean=1;
stats_dozscore = 0;
do_zscore = 0;
showError = 1;

% Time series
    plot_SFC_spectra = 0;
        Ntrials_plot = 0;
    plot_RT = 0;
    
% Stats
    plot_SFC_spectra_stats = 0;
    
% Permutation test results
    plot_Sch_sens = 1;
    
% Scatter plots
    phi_analysis = 0;   % Mainly used in conjunction with figs_batch right now
    plot_regress_alpha_vs_beta = 0;
    
% For dealing with phase plots
    ispaired = 0;
    if plotmode == 4
        %ispaired = 1;
    end
    


    %% Groups Definition
% Setup groups
    
% N_criteria = 3;
% gr = Grp;
% gr.criteria = [2*ones(1,N_criteria)]; gr.criteria_alt = [2*ones(1,N_criteria)]; gr.ctgs = 1;
% i=0; clear group
% i=i+1; group(i) = gr; group(i).ctgs = [1]; group(i).legend = 'All Fast RT';
% i=i+1; group(i) = gr; group(i).ctgs = [2]; group(i).legend = 'All Slow RT';
% 
% i=0; clear group
% % i=i+1; group(i) = gr; group(i).ctgs = [1]; group(i).legend = 'Congruent Fast RT';
% % i=i+1; group(i) = gr; group(i).ctgs = [2]; group(i).legend = 'Congruent Slow RT';
% i=i+1; group(i) = gr; group(i).ctgs = [3]; group(i).legend = 'Incongruent Fast RT';
% i=i+1; group(i) = gr; group(i).ctgs = [4]; group(i).legend = 'Incongruent Slow RT';
% 
% i=0; clear group
% i=i+1; group(i) = gr; group(i).ctgs = [3]; group(i).legend = 'All switch';
% i=i+1; group(i) = gr; group(i).ctgs = [6]; group(i).legend = 'All non-switch';


N_criteria = 3;
gr2 = Grp;
gr2.criteria = [2*ones(1,N_criteria)]; gr2.criteria_alt = [2*ones(1,N_criteria)]; gr2.ctgs = 1;

% Default grouping, SchA vs B
clear group
i=0;
i=i+1; group(i) = gr2; group(i).ctgs = [9]; group(i).legend = 'SchA';
i=i+1; group(i) = gr2; group(i).ctgs = [10]; group(i).legend = 'SchB';            

% Override group special modes (boundary trials, reaction times, switch trials)
group_override = get_group_overrides(fname_suffix_sfc, N_criteria, gr_submode);
if ~isempty(group_override); group = group_override; end


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

% Load the pairs data (coherence)
% SFC
buff_fieldname = ['sfc_' mode2modename(sfc_mode)];
buff_func = @() func_recalc_pairs(stage_range,file_range,sfc_mode);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_sfc,reload_data);
vars_pull(wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc)));
mypairs1 = wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc)).mypairs(:,:,end);
vars_pull(wrkspc_buffer.(buff_fieldname).pairsdat.(stagename(curr_stage_sfc)));
s = wrkspc_buffer.(buff_fieldname).pairsdat.(stagename(curr_stage_sfc)); Ncells_shift1 = s.Ncells_shift; Nelects_shift1 = s.Nelects_shift; clear s

% Percent trials switch (used for sens, but it's really not a useful predictor)
buff_fieldname = ['switchtr'];
buff_func = @() calc_percent_switch(md);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_sfc,reload_data);
vars_pull(wrkspc_buffer.(buff_fieldname));


% Load permutation tests
buff_mode = perm_mode;
buff_fieldname = ['per_' mode2modename(buff_mode)];
buff_func = @() func_recalc_perms(buff_mode,stage_range,file_range);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_permute,reload_data);
% vars_pull(wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc)));
mypairs2 = wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc)).mypairs(:,:,end);
s = wrkspc_buffer.(buff_fieldname).pairsdat.(stagename(curr_stage_sfc)); Ncells_shift2 = s.Ncells_shift; Nelects_shift2 = s.Nelects_shift; clear s

if curr_stage_sfc == 2 || curr_stage_sfc == 3
    % Load bootstrap tests
    buff_mode = boot_mode;
    buff_fieldname = ['boot_' mode2modename(buff_mode)];
    buff_func = @() func_recalc_boot(buff_mode,stage_range,file_range);
    wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_bootstrap,reload_data);
    % vars_pull(wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc)));
    mypairs3 = wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc)).mypairs(:,:,end);
    s = wrkspc_buffer.(buff_fieldname).pairsdat.(stagename(curr_stage_sfc)); Ncells_shift3 = s.Ncells_shift; Nelects_shift3 = s.Nelects_shift; clear s

    if ~test_vars_identical(mypairs1,mypairs2,mypairs3); warning('Mypairs mismatch'); end
    if ~test_vars_identical(Ncells_shift1,Ncells_shift2,Ncells_shift3); warning('Ncells_shift mismatch'); end
    if ~test_vars_identical(Nelects_shift1,Nelects_shift2,Nelects_shift3); warning('Nelects_shift mismatch'); end

end

mypairs = mypairs(:,:,end)';
Ncells_shift = Ncells_shift(:);
Nelects_shift = Nelects_shift(:);
lumap=lumap(:);

% Setup some metadata variables
N_ctgs_base = N_ctgs_base(1);
N_ctgs_extras = N_ctgs_extras(1);
N_ctgs = N_ctgs_base + N_ctgs_extras+1;

% Bad files
[bad_clip, bad_60,fnames_clip] = loadall_bad_data_lfp(curr_stage_badclipping,file_range,md);

