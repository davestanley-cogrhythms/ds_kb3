

% % For loop for saving figs
if ~exist('currfname'); currfname = 'figs_batch_2_0_ctgs'; end
if ~exist('currfigname'); currfigname = 'Fig6'; end
savenames={'fig1','fig2','fig3','fig4','fig5','fig6','fig7','fig8','fig9','fig10','fig11','fig12','fig13','fig14','fig15','fig16','fig17','fig18'};
mydate = datestr(datenum(date),'yy/mm/dd'); mydate = strrep(mydate,'/','');
c=clock;
sp = ['d' mydate '_t' num2str(c(4),'%10.2d') '' num2str(c(5),'%10.2d') '' num2str(round(c(6)),'%10.2d')];
sp = [sp '__' currfname '_' currfigname];
mkdir(fullfile('.',sp));
for i=[8,9,1,2,3,4,5,6]
    figure(i); %ylim([0.75 2])
    title('');
    set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',fullfile('.',sp,savenames{i}))
end
% 
% 
% % % % % % List of frequency peaks (estimated from RT plots of PSD and FFC)
% PSD
% L
% alpha sample NA
% alpha delay 10.74
% beta sample 19.53
% beta delay 17.58
% 
% O
% alpha sample 11.72
% alpha delay 10.74
% beta sample NA
% beta delay NA
% 
% 
% FFC
% L
% alpha sample na
% alpha delay 10.74
% beta sample 17.58
% beta delay 16.6
% 
% O
% alpha sample 
% alpha delay 10.74
% beta sample 17.58
% beta delay 16.6

%% All figures

figs_batch_defaults
addpath('figs_batch_2_0_supporting');


%% Preload all SFC data

% Setup variables
sfc_mode = [22.4413101, 22.4513101, 22.4914103];
perm_mode = [22.4013111, 22.4014111, 22.4413111, 22.4414111, 22.4513111, 22.4714111];
path_matfiles_store = getpath('path_matfiles_store');

% Create wrkspc_buffer structure if not already present
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end

% Load the pairs data (coherence)
for i = 1:length(sfc_mode)
    % SFC
    buff_fieldname = ['sfc_' mode2modename(sfc_mode(i))];
    buff_func = @() func_recalc_pairs(stage_range,file_range,sfc_mode(i));
    wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,0,0);
end

for i = 1:length(perm_mode)
    % Load perm_mode tests
    buff_mode = perm_mode(i);
    buff_fieldname = ['per_' mode2modename(buff_mode)];
    buff_func = @() func_recalc_perms(buff_mode,[2,3],1:79);
    wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,0,0);

end


clear sfc_mode perm_mode


%% Fig 6 - Ctgs1v2 - All - PermuteSpect

clearvars -except wrkspc_buffer fv fv3
vars_pull(fv);
currfigname = 'Fig6';
tic
opts_PM3Dcs.do_sgolay = 0;
opts_PM3Dcs.stats_mode = 2;
opts_PM3Dcs.remove_dependent = 1;

do_both_animals=1;
ffc_sfc_psd_mode = 1;
switch ffc_sfc_psd_mode
    case 1
        %sfc_mode=22.4013101;    
        perm_mode=22.4014111;
        script_name = 'plots_ffc_vs_Sch';
        freqband_stats = freqband_stats_lbeta;
        force_freqband_stats = 0;
    case 2
        %sfc_mode=2.2013100;     
        perm_mode=2.2015111;
        script_name = 'plots_preferred';
        freqband_stats = freqband_stats_hbeta;
        force_freqband_stats = 0;
    case 3
        %sfc_mode=41.6013101;
        perm_mode=41.6014111;       
        script_name = 'plots_CFC';
        freqband_stats = freqband_stats_alpha;
        force_freqband_stats = 1;
end

% [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    
%exclude_clipping = 1;  warning('Including clip!');
exclude_60 = 0;

permdat2pls = 0;            % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    do_diff = 0;            % Actually, don't normalize by null distribution
    sort_pls = 0;
perm2pls = 1;               % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
    perm2pls_dophi = 0;
    
plot_SFC_spectra=0;

remove_dependent = 0;

if force_freqband_stats;
    if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end;
end
% Choose animals & stage
excludeL = 0;
if do_both_animals; excludeO = 0; warning('Doing all animals'); else excludeO=1; end
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
groupls = out.group;
plsls = out.pls_stats;
bad_ls = out.bad_any;
clear out

if force_freqband_stats;
    if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = freqband_stats_alpha; end;
end
% Choose animals & stage
excludeL = 0; 
if do_both_animals; excludeO = 0; warning('Doing all animals'); else excludeO=1; end
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
groupld = out.group;
plsld = out.pls_stats;
bad_ld = out.bad_any;
clear out


if ~do_both_animals
    if force_freqband_stats;
        if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = freqband_stats_lbeta; end;
    end
    % Choose animals & stage
    excludeL = 1; 
    excludeO = 0; 
    curr_stage_sfc=2;
        curr_stage_badclipping = curr_stage_sfc;
    push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
    groupos = out.group;
    plsos = out.pls_stats;
    bad_os = out.bad_any;
    clear out

    if force_freqband_stats;
        if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = freqband_stats_alpha; end;
    end
    % Choose animals & stage
    excludeL = 1; 
    excludeO = 0; 
    curr_stage_sfc=3;
        curr_stage_badclipping = curr_stage_sfc;
    push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
    groupod = out.group;
    plsod = out.pls_stats;
    bad_od = out.bad_any;
    clear out
end
toc

%
% inds = [1:length(groupls)];
% inds = [1,2,3,4,9,10];
% names = {'SchA Pref','SchA NonPref','SchB Pref','SchB NonPref','Pref Sch','NonPref Sch'};

inds = [1,2,5];
% names = {'SchA Pref - NonPref','SchB Pref - NonPref','Pref Sch - NonPref Sch'};
names = {'Cat vs Dog','Fat vs Thin','SchA vs B'};

% inds = [1:5];
% names = {'Ctg1 vs 2','Ctg3 vs 4','Ctg1 vs 2 irr','Ctg3 vs 4 irr','SchA vs B'};

% inds = [1:3];
% names = {'Ctg 1 vs 2','Ctg 3 vs 4','Sch A vs B'};
% names = {'SchA Pref - NonPref','SchB Pref - NonPref','Pref Sch - NonPref Sch'};

[groupls(inds).legend] = names{:};
[groupld(inds).legend] = names{:};
[groupos(inds).legend] = names{:};
[groupod(inds).legend] = names{:};

%
i=0;
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupls(inds) ],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupld(inds) ],[],opts_PM3Dcs,[]);
if ~do_both_animals
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupos(inds) ],[],opts_PM3Dcs,[]);
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupod(inds) ],[],opts_PM3Dcs,[]);
end

%%
% opts_PSC.hmask = blkdiag(true(2),true(2),true(2),true(2),true(2),true(2),true(2),true(2),true(2),true(2),true(2),true(2));
% i=i+1; figure; [h1, h, p, out.PSC] = plot_stats_custstruct([groupls(inds) groupld(inds) groupos(inds) groupod(inds)],opts_PSC);

opts_PSC.hmask = blkdiag(true(3),true(3),true(3),true(3));
i=i+1; figure; [h1, h, p, out.PSC] = plot_stats_custstruct([groupls(inds) groupld(inds) groupos(inds) groupod(inds)],opts_PSC);


% Save figures
% for i=1:4; figure(i); set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i}); end

%% Monkey L Sample
% General variables for all cases
opts_PPS_temp = opts_PPS;
opts_PPS_temp.hide_all_plots = 1;
opts_PPS_temp.do_bh0=0;
opts_PPS_temp.do_phi=perm2pls_dophi;

plot_scatter_both_animals = 1;

bad_curr = bad_ls;
pls_curr = plsls; 

sp_ls = [];
for i = 1:length(groupls)
    curr_stage_sfc = 2;
    opts_PPS_temp.sens_chosen_condition = i;
    
    
    [PPS_temp] = plot_permute_scatter(...
        wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
        [], ...
        bad_curr,freqband_stats, opts_PPS_temp);
    sp_ls = [sp_ls PPS_temp.q4_perm(:)];
    clear PPS_temp
end

warning('Might need to still fix this');
N_criteria = 5;
gr_template = Grp;
gr_template.criteria = [2*ones(1,N_criteria)]; gr_template.criteria_alt = []; gr_template.criteria_sfc = []; gr_template.ctgs = 1;
gr_template.xlims_desired = [0 120];


clear group
i=0;
% Ctg A Sensitive - Pref vs non-preferred
i=i+1; group(i) = gr_template; group(i).criteria=[1 0 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch A';
i=i+1; group(i) = gr_template; group(i).criteria=[0 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch B';
i=i+1; group(i) = gr_template; group(i).criteria=[1 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Both';

for i = 1:length(group); [group(i).cells, group(i).numcells]= get_grouped_cells(group(i),[sp_ls]); end

if plot_scatter_both_animals
    figure; plot(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),'kx','LineWidth',2); hold on; plot_scattergroups(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),group,~bad_curr); xlabel('FFC Ctg1-Ctg2'); ylabel('FFC Ctg3-Ctg4')
else
    figure; plot(pls_curr(:,1),pls_curr(:,2),'k.'); hold on; plot_scattergroups(pls_curr(:,1),pls_curr(:,2),group); xlabel('FFC Ctg1-Ctg2'); ylabel('FFC Ctg3-Ctg4')
end
title('Monkey L Sample');

% Monkey L Delay
bad_curr = bad_ld;
pls_curr = plsld; 

sp_ls = [];
for i = 1:length(groupls)
    curr_stage_sfc = 3;
    opts_PPS_temp.sens_chosen_condition = i;
    
    
    [PPS_temp] = plot_permute_scatter(...
        wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
        [], ...
        bad_curr,freqband_stats, opts_PPS_temp);
    sp_ls = [sp_ls PPS_temp.q4_perm(:)];
    clear PPS_temp
end

N_criteria = 5;
gr_template = Grp;
gr_template.criteria = [2*ones(1,N_criteria)]; gr_template.criteria_alt = []; gr_template.criteria_sfc = []; gr_template.ctgs = 1;
gr_template.xlims_desired = [0 120];


clear group
i=0;
% Ctg A Sensitive - Pref vs non-preferred
i=i+1; group(i) = gr_template; group(i).criteria=[1 0 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch A';
i=i+1; group(i) = gr_template; group(i).criteria=[0 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch B';
i=i+1; group(i) = gr_template; group(i).criteria=[1 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Both';

for i = 1:length(group); [group(i).cells, group(i).numcells]= get_grouped_cells(group(i),[sp_ls]); end


if plot_scatter_both_animals
    figure; plot(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),'kx','LineWidth',2); hold on; plot_scattergroups(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),group,~bad_curr); xlabel('FFC Ctg1-Ctg2'); ylabel('FFC Ctg3-Ctg4')
else
    figure; plot(pls_curr(:,1),pls_curr(:,2),'k.'); hold on; plot_scattergroups(pls_curr(:,1),pls_curr(:,2),group); xlabel('FFC Ctg1-Ctg2'); ylabel('FFC Ctg3-Ctg4')
end
title('Monkey L Delay');

%% Monkey O Sample
bad_curr = bad_os;
pls_curr = plsos; 

sp_ls = [];
for i = 1:length(groupls)
    curr_stage_sfc = 2;
    opts_PPS_temp.sens_chosen_condition = i;
    
    
    [PPS_temp] = plot_permute_scatter(...
        wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
        [], ...
        bad_curr,freqband_stats, opts_PPS_temp);
    sp_ls = [sp_ls PPS_temp.q4_perm(:)];
    clear PPS_temp
end

N_criteria = 5;
gr_template = Grp;
gr_template.criteria = [2*ones(1,N_criteria)]; gr_template.criteria_alt = []; gr_template.criteria_sfc = []; gr_template.ctgs = 1;
gr_template.xlims_desired = [0 120];


clear group
i=0;
% Ctg A Sensitive - Pref vs non-preferred
i=i+1; group(i) = gr_template; group(i).criteria=[1 0 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch A';
i=i+1; group(i) = gr_template; group(i).criteria=[0 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch B';
i=i+1; group(i) = gr_template; group(i).criteria=[1 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Both';

for i = 1:length(group); [group(i).cells, group(i).numcells]= get_grouped_cells(group(i),[sp_ls]); end


if plot_scatter_both_animals
    figure; plot(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),'kx','LineWidth',2); hold on; plot_scattergroups(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),group,~bad_curr); xlabel('FFC Ctg1-Ctg2'); ylabel('FFC Ctg3-Ctg4')
else
    figure; plot(pls_curr(:,1),pls_curr(:,2),'k.'); hold on; plot_scattergroups(pls_curr(:,1),pls_curr(:,2),group); xlabel('FFC Ctg1-Ctg2'); ylabel('FFC Ctg3-Ctg4')
end
title('Monkey O Sample');

% Monkey O Delay
bad_curr = bad_od;
pls_curr = plsod; 

sp_ls = [];
for i = 1:length(groupls)
    curr_stage_sfc = 3;
    opts_PPS_temp.sens_chosen_condition = i;
    
    
    [PPS_temp] = plot_permute_scatter(...
        wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
        [], ...
        bad_curr,freqband_stats, opts_PPS_temp);
    sp_ls = [sp_ls PPS_temp.q4_perm(:)];
    clear PPS_temp
end

N_criteria = 5;
gr_template = Grp;
gr_template.criteria = [2*ones(1,N_criteria)]; gr_template.criteria_alt = []; gr_template.criteria_sfc = []; gr_template.ctgs = 1;
gr_template.xlims_desired = [0 120];


clear group
i=0;
% Ctg A Sensitive - Pref vs non-preferred
i=i+1; group(i) = gr_template; group(i).criteria=[1 0 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch A';
i=i+1; group(i) = gr_template; group(i).criteria=[0 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch B';
i=i+1; group(i) = gr_template; group(i).criteria=[1 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Both';

for i = 1:length(group); [group(i).cells, group(i).numcells]= get_grouped_cells(group(i),[sp_ls]); end

if plot_scatter_both_animals
    figure; plot(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),'kx','LineWidth',2); hold on; plot_scattergroups(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),group,~bad_curr); xlabel('FFC Ctg1-Ctg2'); ylabel('FFC Ctg3-Ctg4')
else
    figure; plot(pls_curr(:,1),pls_curr(:,2),'k.'); hold on; plot_scattergroups(pls_curr(:,1),pls_curr(:,2),group); xlabel('FFC Ctg1-Ctg2'); ylabel('FFC Ctg3-Ctg4')
end
title('Monkey O Delay');

% Fig 7 - Ctgs Scatterplots


%% Monkey O alpha vs beta scatter

clear group
freqband_stats = freqband_stats_alpha;
% Choose animals & stage
excludeL = 1; 
excludeO = 0; 
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
plsoda = out.pls_stats;
clear out


bad_curr = bad_od;
pls_curr = plsod; 

sp_ls = [];
for i = 1:length(groupls)
    curr_stage_sfc = 3;
    opts_PPS_temp.sens_chosen_condition = i;
    
    
    [PPS_temp] = plot_permute_scatter(...
        wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
        [], ...
        bad_curr,freqband_stats, opts_PPS_temp);
    sp_ls = [sp_ls PPS_temp.q4_perm(:)];
    clear PPS_temp
end

N_criteria = 5;
gr_template = Grp;
gr_template.criteria = [2*ones(1,N_criteria)]; gr_template.criteria_alt = []; gr_template.criteria_sfc = []; gr_template.ctgs = 1;
gr_template.xlims_desired = [0 120];


clear group
i=0;
% Ctg A Sensitive - Pref vs non-preferred
i=i+1; group(i) = gr_template; group(i).criteria=[1 0 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch A';
i=i+1; group(i) = gr_template; group(i).criteria=[0 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch B';
i=i+1; group(i) = gr_template; group(i).criteria=[1 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Both';

for i = 1:length(group); [group(i).cells, group(i).numcells]= get_grouped_cells(group(i),[sp_ls]); end

figure; plot(pls_curr(~bad_curr,1),plsoda(~bad_curr,1),'kx','LineWidth',2); hold on; plot_scattergroups(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),group,~bad_curr); xlabel('Beta FFC Ctg1-Ctg2'); ylabel('FFC Ctg1-Ctg2')
figure; plot(pls_curr(~bad_curr,1),plsoda(~bad_curr,2),'kx','LineWidth',2); hold on; plot_scattergroups(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),group,~bad_curr); xlabel('Beta FFC Ctg1-Ctg2'); ylabel('FFC Ctg3-Ctg4')
figure; plot(pls_curr(~bad_curr,2),plsoda(~bad_curr,1),'kx','LineWidth',2); hold on; plot_scattergroups(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),group,~bad_curr); xlabel('Beta FFC Ctg3-Ctg4'); ylabel('FFC Ctg1-Ctg2')
figure; plot(pls_curr(~bad_curr,2),plsoda(~bad_curr,2),'kx','LineWidth',2); hold on; plot_scattergroups(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),group,~bad_curr); xlabel('Beta FFC Ctg3-Ctg4'); ylabel('FFC Ctg3-Ctg4')


title('Monkey O Delay');


%% Fig 6b - Ctgs Rel vs Irrel
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv);
currfigname = 'Fig6b';

remove_dependent = 1;

ffc_sfc_psd_mode = 1;
switch ffc_sfc_psd_mode
    case 1
        %sfc_mode=22.4013101;    
        perm_mode=22.4014111;           % NOTE : HAVE TO USE THIS ONE - SINCE SOME DATA MIGHT BE MISSING.
        perm_mode2=22.4014111;          % DITO.
        script_name = 'plots_ffc_vs_Sch';
        freqband_stats = freqband_stats_lbeta;
        force_freqband_stats = 0;
    case 2
        %sfc_mode=2.2014101;     
        perm_mode=2.2015111;
        perm_mode2=2.2015111;
        script_name = 'plots_preferred';
        freqband_stats = freqband_stats_hbeta;
        force_freqband_stats = 0;
    case 3
        %sfc_mode=41.6014101;
        perm_mode=41.6014111;       
        perm_mode2=41.6014111;       
        script_name = 'plots_CFC';
        freqband_stats = freqband_stats_alpha;
        force_freqband_stats = 1;
end

% Load perm_mode2 tests
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end
path_matfiles_store = getpath('path_matfiles_store');
buff_mode = perm_mode2;
buff_fieldname = ['per_' mode2modename(buff_mode)];
buff_func = @() func_recalc_perms(buff_mode,[2,3],1:79);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,0,0);



% [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
%%
%exclude_clipping = 1;  warning('Including clip!');
exclude_60 = 0;

permdat2pls = 1;            % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    do_diff = 0;            % Actually, don't normalize by null distribution
    sort_pls = 0;
perm2pls = 0;               % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
    perm2pls_dophi = 0;
    
plot_SFC_spectra=0;



%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 0;
% Choose stage
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
%
Fig_6b_Ctgs_Rel_Irrel(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)


%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 0;
% Choose stage
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
%
Fig_6b_Ctgs_Rel_Irrel(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)

%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 1;
% Choose stage
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
%
Fig_6b_Ctgs_Rel_Irrel(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)


%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 1;
% Choose stage
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
%
Fig_6b_Ctgs_Rel_Irrel(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)


%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 1;
excludeO = 0;
% Choose stage
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
%
Fig_6b_Ctgs_Rel_Irrel(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)


%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 1;
excludeO = 0;
% Choose stage
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
%
Fig_6b_Ctgs_Rel_Irrel(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)


%% Fig 7 - Ctg Selectivity - Boundary vs Nonboundary
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv);
currfigname = 'Fig7';

do_both_animals=0;
ffc_sfc_psd_mode = 1;
switch ffc_sfc_psd_mode
    case 1
        %sfc_mode=22.4013101;    
        perm_mode=22.4714111;
        perm_mode2=22.4014111;              % Used to determine sorting order of preferred / non-preferred
        script_name = 'plots_ffc_vs_Sch';
        freqband_stats = freqband_stats_lbeta;
        force_freqband_stats = 0;
    case 2
        %sfc_mode=2.2013100;     
        perm_mode=2.2715111;
        perm_mode2=2.2015111;
        script_name = 'plots_preferred';
        freqband_stats = freqband_stats_hbeta;
        force_freqband_stats = 0;
    case 3
        %sfc_mode=41.6013101;
        perm_mode=41.6714111;       
        perm_mode2=41.6014111; 
        script_name = 'plots_CFC';
        freqband_stats = freqband_stats_alpha;
        force_freqband_stats = 1;
end

% Load perm_mode2 tests
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end
path_matfiles_store = getpath('path_matfiles_store');
buff_mode = perm_mode2;
buff_fieldname = ['per_' mode2modename(buff_mode)];
buff_func = @() func_recalc_perms(buff_mode,[2,3],1:79);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,0,0);


% [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    
%exclude_clipping = 1;  warning('Including clip!');
exclude_60 = 0;

permdat2pls = 1;            % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    do_diff = 0;            % Actually, don't normalize by null distribution
    sort_pls = 0;
perm2pls = 0;               % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
    perm2pls_dophi = 0;
    
plot_SFC_spectra=0;


% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
if do_both_animals;
    excludeO = 0; warning('Doing all animals'); else excludeO=1;
end
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
Fig_7_Ctg_Selectivity(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)

%%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
if do_both_animals;
    excludeO = 0; warning('Doing all animals'); else excludeO=1;
end
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
Fig_7_Ctg_Selectivity(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)



%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 1;
if do_both_animals;
    excludeO = 0; warning('Doing all animals'); else excludeO=0;
end
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
Fig_7_Ctg_Selectivity(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)

%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 1;
if do_both_animals;
    excludeO = 0; warning('Doing all animals'); else excludeO=0;
end
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
Fig_7_Ctg_Selectivity(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)


%% Fig 8a - Do ctg specialists also have higher bndry response? (Ctgsetli mode 0 vs 4)
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv);
currfigname = 'Fig8a';

remove_dependent = 0;
ffc_sfc_psd_mode = 1;
switch ffc_sfc_psd_mode
    case 1
        sfc_mode=22.4413101;    
        perm_mode=22.4413111;
        perm_mode2=22.4014111;              % Used to determine sorting order of preferred / non-preferred
        script_name = 'plots_ffc_vs_Sch';
        freqband_stats = freqband_stats_lbeta;
        force_freqband_stats = 0;
    case 2
        sfc_mode=2.2413100;     
        perm_mode=2.2415111;
        perm_mode2=2.2015111;
        script_name = 'plots_preferred';
        freqband_stats = freqband_stats_hbeta;
        force_freqband_stats = 0;
    case 3
        sfc_mode=41.6413101;
        perm_mode=41.6413111;       
        perm_mode2=41.6014111; 
        script_name = 'plots_CFC';
        freqband_stats = freqband_stats_alpha;
        force_freqband_stats = 1;
end


% Load perm_mode2 tests
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end
path_matfiles_store = getpath('path_matfiles_store');
buff_mode = perm_mode2;
buff_fieldname = ['per_' mode2modename(buff_mode)];
buff_func = @() func_recalc_perms(buff_mode,[2,3],1:79);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,0,0);



% [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    
%exclude_clipping = 1;  warning('Including clip!');
exclude_60 = 0;

permdat2pls = 0;            % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    do_diff = 0;            % Actually, don't normalize by null distribution
    sort_pls = 0;
perm2pls = 0;               % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
    perm2pls_dophi = 0;
    
plot_SFC_spectra=0;


%%
% Sample/Sample
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 0;
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
mygroup = Fig_8a_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi);
mygroupbs = mygroup;


% Delay/Delay
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 0;
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
mygroup = Fig_8a_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi);
mygroupbd = mygroup;

%%
% Sample/Delay
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 0;
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
curr_stage_sfc=3;
mygroup = Fig_8a_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi);
mygroupbs = mygroup;


%%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 1;
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
mygroup = Fig_8a_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi);
mygroupls = mygroup;

%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 1;
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
mygroup = Fig_8a_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi);
mygroupld = mygroup;

%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 1;
excludeO = 0;

curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
mygroup = Fig_8a_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi);
mygroupos = mygroup;

%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 1;
excludeO = 0;

curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
mygroup = Fig_8a_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi);
mygroupod = mygroup;

% Plotting
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroupbs],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroupbd],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroupls],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroupld],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroupos],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroupod],[],opts_PM3Dcs,[]);

% Stats
opts_PSC.hmask = blkdiag(true(2),true(2));
i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroupbs mygroupbd],opts_PSC);

opts_PSC.hmask = blkdiag(true(2),true(2),true(2),true(2));
i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroupls mygroupld mygroupos mygroupod],opts_PSC);


%% Fig 8b - Ctgsetli mode 0 vs mode 5
% Has 3 modes:
% 1) Do specialists have higher boundary response than non-specialists? (Fig8b_mode=1)
% 2) Also, what is specialist boundary rel vs irrel? (Fig8b_mode=2)
% 3) Scatterplot of bA-NbA vs bB-NbB with specialist groupings  (Fig8b_mode=3)
% 4) Scatterplot of Ctg sensitivity vs boundary sensitivity. (Fig8b_mode=4)
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv);
currfigname = 'Fig8b';
rel_vs_irrel = 0;
Fig8b_mode = 4;

do_both_animals=0;
ffc_sfc_psd_mode = 1;
switch ffc_sfc_psd_mode
    case 1
        sfc_mode=22.4513101;    
        perm_mode=22.4513111;
        perm_mode2=22.4014111;              % Used to determine sorting order of preferred / non-preferred
        script_name = 'plots_ffc_vs_Sch';
        freqband_stats = freqband_stats_lbeta;
        force_freqband_stats = 0;
    case 2
        sfc_mode=2.2513100;     
        perm_mode=2.2515111;
        perm_mode2=2.2015111;
        script_name = 'plots_preferred';
        freqband_stats = freqband_stats_hbeta;
        force_freqband_stats = 0;
    case 3
        sfc_mode=41.6513101;
        perm_mode=41.6513111;       
        perm_mode2=41.6014111; 
        script_name = 'plots_CFC';
        freqband_stats = freqband_stats_alpha;
        force_freqband_stats = 1;
end


% Load perm_mode2 tests
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end
path_matfiles_store = getpath('path_matfiles_store');
buff_mode = perm_mode2;
buff_fieldname = ['per_' mode2modename(buff_mode)];
buff_func = @() func_recalc_perms(buff_mode,[2,3],1:79);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,0,0);



% [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    
%exclude_clipping = 1;  warning('Including clip!');
exclude_60 = 0;

permdat2pls = 0;            % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    do_diff = 0;            % Actually, don't normalize by null distribution
    sort_pls = 0;
perm2pls = 0;               % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
    perm2pls_dophi = 0;
    
plot_SFC_spectra=0;

%%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
if do_both_animals;
    excludeO = 0; warning('Doing all animals'); else excludeO=1;
end
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
Fig_8b_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi,Fig8b_mode);

%%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
if do_both_animals;
    excludeO = 0; warning('Doing all animals'); else excludeO=1;
end
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
Fig_8b_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi,Fig8b_mode);



%%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 1;
excludeO = 0;
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
Fig_8b_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi,Fig8b_mode);

%%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 1;
excludeO = 0;
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
Fig_8b_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi,Fig8b_mode);



%% Fig 8c - Do specialist cells exhibit enhanced boundary response for their preferred boundary?
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv);
currfigname = 'Fig8c';

ffc_sfc_psd_mode = 1;
switch ffc_sfc_psd_mode
    case 1
        sfc_mode=22.4514101;    
        perm_mode=22.4513111;
        perm_mode2=22.4014111;              % Used to determine sorting order of preferred / non-preferred
        script_name = 'plots_ffc_vs_Sch';
        freqband_stats = freqband_stats_lbeta;
        force_freqband_stats = 0;
    case 2
        sfc_mode=2.2513100;     
        perm_mode=2.2515111;
        perm_mode2=2.2015111;
        script_name = 'plots_preferred';
        freqband_stats = freqband_stats_hbeta;
        force_freqband_stats = 0;
    case 3
        sfc_mode=41.6513101;
        perm_mode=41.6513111;       
        perm_mode2=41.6014111; 
        script_name = 'plots_CFC';
        freqband_stats = freqband_stats_alpha;
        force_freqband_stats = 1;
end


% Load perm_mode2 tests
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end
path_matfiles_store = getpath('path_matfiles_store');
buff_mode = perm_mode2;
buff_fieldname = ['per_' mode2modename(buff_mode)];
buff_func = @() func_recalc_perms(buff_mode,[2,3],1:79);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,0,0);



% [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    
%exclude_clipping = 1;  warning('Including clip!');
exclude_60 = 0;

permdat2pls = 0;            % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    do_diff = 0;            % Actually, don't normalize by null distribution
    sort_pls = 0;
perm2pls = 0;               % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
    perm2pls_dophi = 0;
    
plot_SFC_spectra=0;
remove_dependent=0;

%%

% Sample/Sample
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 0;
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
Fig_8c_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)
%%
% Sample/Delay
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 0;
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
curr_stage_sfc=2;
Fig_8c_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)



%%
% Delay/Delay
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 0;
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
Fig_8c_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)


%%

% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 1;
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
Fig_8c_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)

%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 1;
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
Fig_8c_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)



%%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 1;
excludeO = 0;
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
Fig_8c_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)

%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 1;
excludeO = 0;
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
Fig_8c_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)




%% Fig 9 - FFC Royfig
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv);
currfigname = 'Fig9';

ffc_sfc_psd_mode = 1;
switch ffc_sfc_psd_mode
    case 1
        sfc_mode=22.4914103;    
        perm_mode=22.4014111;
        perm_mode2=22.4513111;
        script_name = 'plots_ffc_vs_Sch';
        freqband_stats = freqband_stats_lbeta;
        force_freqband_stats = 0;
    case 2
        sfc_mode=2.2914103;     
        perm_mode=2.2015111;
        perm_mode2=2.2615111;
        script_name = 'plots_preferred';
        freqband_stats = freqband_stats_hbeta;
        force_freqband_stats = 0;
    case 3
        sfc_mode=41.6914103;
        perm_mode=41.6014111;       
        perm_mode2=41.6614111;       
        script_name = 'plots_CFC';
        freqband_stats = freqband_stats_alpha;
        force_freqband_stats = 1;
end

% Load perm_mode2 tests
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end
path_matfiles_store = getpath('path_matfiles_store');
buff_mode = perm_mode2;
buff_fieldname = ['per_' mode2modename(buff_mode)];
buff_func = @() func_recalc_perms(buff_mode,[2,3],1:79);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,0,0);


% [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);

%exclude_clipping = 1;  warning('Including clip!');
exclude_60 = 0;

permdat2pls = 0;            % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    do_diff = 0;            % Actually, don't normalize by null distribution
    sort_pls = 0;
perm2pls = 0;               % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
    perm2pls_dophi = 0;
    
plot_SFC_spectra=0;


% Choose stage
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;

%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 0;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
%
Fig_9_Royfig(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)

%%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 0;
excludeO = 1;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
%
Fig_9_Royfig(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)

%%
%
% Choose animals & stage
if force_freqband_stats; if ffc_sfc_psd_mode == 3; warning('Forcing freqband_stats'); freqband_stats = [26 28]; end; end
excludeL = 1;
excludeO = 0;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
Fig_9_Royfig(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)


