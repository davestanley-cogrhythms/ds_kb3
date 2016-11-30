
%% Stage
% clear
clc
run_setdefaultfig
addpath(genpath('./funcs_plots_preferred/'));

% Data loading parameters
recalc_md = 0;
recalc_cfc = 0;
recalc_psd = 0;
recalc_psd_permute = 1;
recalc_pac = 0;
recalc_pac_permute = 0;
recalc_bootstrap = 0;
reload_data = 0;
if get_iscromer; stage_range = [2,3];
else stage_range = [2,3];
end
% sfc_mode=41.6613101;
% perm_mode=41.6613111;
cfc_mode=40.60131;  % (amp-amp coupling)
pac_mode=43.601210;
pac_mode_permute = 43.601311;

% perm_mode = sfc_mode+0.000001;   % Make sure these are aligned

group_mode = 1;
gr_submode = 1;             % Which of the ctgs to plot from the bndry, RT modes (i.e. congruent/incongruent; 50 or 60% morphs

[~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
fprintf(['Chosen PSD mode: ' fname_suffix_sfc ' \n']);
[~, cfc_subgroups] = decode_sfc_mode(cfc_mode);
[fname_suffix_cfc] = build_sfcmode(cfc_mode, cfc_subgroups);


% Load preferred units
    % Stage selection

    curr_stage = 2;
    curr_stage_alt = 3;
    curr_stage_sfc = 2;
    curr_stage_badclipping = curr_stage_sfc;

    % Unit exclusion
    exclude_clipping = 1;
    exclude_60 = 0;
    excludeL = 1;
    excludeO = 0;
    
    
    % Grouping switches
    groups_separate_monkeys = 0;
    preferred_mode = 3;     % 1 for statistics (ttest, not working in plots_preferred_unitpairs.m); 2 for raw threshold; 3 for quantile
    
    % CFC parameters
%     fix these indicies - need to reflect current frequency vector! - otherwise it messes up when switch curr_stage from 2 to 3
    freqband_stats = [1 4];
%     freqband_corrcoef = [11];    % Measure correlations with this frequency
    freqband_stats = [10 12];
%     freqband_stats = [16 20];
    
    plotmode = 1;               % 1-Corrcef; 2-Corrcoef P; 3-S1ave
    perm2pls = 1;               % Instead of showing raw SFC values, show % cells successfully passing permutation test
        perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
        perm2pls_dophi = 0;
    

    opts_PM3Dcs.paperfig_mode=0;
    opts_PM3Dcs.do_sgolay=0;
    opts_PSC.paperfig_mode=0;
    opts_PPS.paperfig_mode=0;


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




%% Plot Switches

% General
domean=1;
stats_dozscore = 0;
do_zscore = 0;
showError = 1;

% Plotting toggles
% Time series
    plot_SFC_spectra = 1;
    plot_PAC = 0;
    plot_PAC_sens = 0;
    
% Stats
    plot_SFC_spectra_stats = 0;


% Scatter plots
    plot_regress_alpha_vs_beta = 0;

% Permutation test results
    plot_Sch_sens = 0;
    
%% Groups Definition
% Setup groups
clear gr2;



switch group_mode
    case 1
        do_gamma = 0;
        N_criteria = 3;
        gr2 = Grp;
        gr2.criteria = [2*ones(1,N_criteria)]; gr2.criteria_alt = [2*ones(1,N_criteria)]; gr2.ctgs = 11;


        clear group
        i=0;
        i=i+1; group(i) = gr2; group(i).criteria = [1 0 0]; group(i).legend = 'Delta';
        i=i+1; group(i) = gr2; group(i).criteria = [0 1 0]; group(i).legend = 'Alpha';
        i=i+1; group(i) = gr2; group(i).criteria = [0 0 1]; group(i).legend = 'Beta';
        % i=i+1; group(i) = gr2; group(i).criteria = [1 1 1]; group(i).legend = 'all';      % No electrodes
        i=i+1; group(i) = gr2; group(i).criteria = [0 0 0]; group(i).legend = 'none';
%         i=i+1; group(i) = gr2; group(i).criteria = [1 1 0]; group(i).legend = 'delta-alpha';
        % i=i+1; group(i) = gr2; group(i).criteria = [1 0 1]; group(i).legend = 'delta-beta'; % No electrodes
%         i=i+1; group(i) = gr2; group(i).criteria = [0 1 1]; group(i).legend = 'beta-alpha';

        
        
    case 2
        do_gamma = 1;
        N_criteria = 4;
        gr2 = Grp;
        gr2.criteria = [2*ones(1,N_criteria)]; gr2.criteria_alt = [2*ones(1,N_criteria)]; gr2.ctgs = 11;


        clear group
        i=0;
        i=i+1; group(i) = gr2; group(i).criteria = [1 0 0 0]; group(i).legend = 'Delta';
        i=i+1; group(i) = gr2; group(i).criteria = [0 1 0 0]; group(i).legend = 'Alpha';
        i=i+1; group(i) = gr2; group(i).criteria = [0 0 1 0]; group(i).legend = 'Beta';
        i=i+1; group(i) = gr2; group(i).criteria = [0 0 0 1]; group(i).legend = 'Gamma';
        i=i+1; group(i) = gr2; group(i).criteria = [0 0 0 0]; group(i).legend = 'none';
end
  

% Override group special modes (boundary trials, reaction times, switch trials)
group_override = get_group_overrides(fname_suffix_sfc, N_criteria, gr_submode);
if ~isempty(group_override); group = group_override; end


if groups_separate_monkeys
    group_M(1:2) = Grp;
    group_M(1).criteria_sfc(2) = 1; group_M(1).legend = 'Monkey1';
    group_M(2).criteria_sfc(2) = 0; group_M(2).legend = 'Monkey2';
    
    %group = group_OUTER(group,group_M); % Separate all groups according to Monkeys
    group = group_M;        % Plot just 2 separate Monkeys
end

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


% Load CFC data (amplitude-amplitude coupling)
% buff_fieldname = ['cfc_' mode2modename(cfc_mode)];
% buff_func = @() func_recalc_mode40 (cfc_mode,stage_range,file_range);
% wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_cfc,reload_data);
% vars_pull(wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc)));
% s = wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc));
% Rtr = s.Cave;
% f = f(1,:);
% clear Cave

% Load PSD data
buff_fieldname_psd = ['psd_' mode2modename(sfc_mode)];
buff_func = @() func_recalc_mode41 (sfc_mode,stage_range,file_range);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname_psd, buff_func,recalc_psd,reload_data);
s = wrkspc_buffer.(buff_fieldname_psd).(stagename(curr_stage_sfc));
S1ave = s.Cave;
f2 = s.f(1,:);
funames1D = s.funames1D;
clear s;

% Load PACdata
buff_fieldname_pac = ['pac_' mode2modename(pac_mode)];
buff_func = @() func_recalc_modePAC (pac_mode,stage_range,file_range);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname_pac, buff_func,recalc_pac,reload_data);
if curr_stage_sfc == 3
    s = wrkspc_buffer.(buff_fieldname_pac).(stagename(curr_stage_sfc));
    Pac = s.Cave;
    pPac = s.p;
    fPac1 = s.f(1,:);
    fPac2 = s.f2(1,:);
    clear s;
end


% Load PSD permutation tests (should work with any case > 40!)
buff_mode = perm_mode;
buff_fieldname = ['per_' mode2modename(buff_mode)];
buff_func = @() func_recalc_perms(buff_mode,stage_range,file_range);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_psd_permute,reload_data);
% vars_pull(wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc)));
mypairs2 = wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc)).mypairs(:,:,end);
s = wrkspc_buffer.(buff_fieldname).pairsdat.(stagename(curr_stage_sfc)); Ncells_shift2 = s.Ncells_shift; Nelects_shift2 = s.Nelects_shift; clear s



% Load CFC/PAC permutation tests (should work with any case > 40!)
buff_fieldname = ['per_' mode2modename(pac_mode_permute)];
buff_func = @() func_recalc_perms_CFC(pac_mode_permute,stage_range,file_range);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_pac_permute,reload_data);
if curr_stage_sfc == 3
    s = wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc));
    dPac9 = s.Cave1;
    dPac10 = s.Cave2;
    pdPac = s.pvals_Cave;
    fdPac1 = s.f(1,:);
    fdPac2 = s.f2(1,:);
    clear s;
end



% Load specialization indices
load(fullfile(getpath('path_buffer_specialization_stats'),'specialization.mat'));
statsnames = stats.(stagename(curr_stage_sfc)).stats(1).m.ctgsetnames;
statsFR = stats.(stagename(curr_stage_sfc)).mu_arr;
statsFR = squeeze(statsFR);


% Setup some metadata variables
% N_ctgs_base = N_ctgs_base(1);
% N_ctgs_extras = N_ctgs_extras(1);
% N_ctgs = N_ctgs_base + N_ctgs_extras+1;

% Bad files
[bad_clip, bad_60,fnames_clip] = loadall_bad_data_lfp(curr_stage_badclipping,file_range,md);
