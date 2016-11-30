
% deleteme

% % % % Random leftover things
% 
% % Old save command
% saveas(gcf,savenames{i},savetype);
% 
% % For loop for saving figs
savenames={'fig1','fig2','fig3','fig4','fig5','fig6','fig7','fig8','fig9','fig10','fig11','fig12','fig13','fig14','fig15','fig16','fig17','fig18'};
mydate = datestr(datenum(date),'yy/mm/dd'); mydate = strrep(mydate,'/','');
c=clock;
savepath = ['fig_d' mydate '_t' num2str(c(4),'%10.2d') '' num2str(c(5),'%10.2d') '' num2str(round(c(6)),'%10.2d')];
mkdir(fullfile('.',savepath));
for i=[1:2]
    figure(i); %ylim([0 15])
    set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',fullfile('.',savepath,savenames{i}))
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


%% Fig 6 - Ctgs1v2 - All - PermuteSpect
clearvars -except wrkspc_buffer fv
vars_pull(fv);

opts_PM3Dcs.do_sgolay = 0;
opts_PM3Dcs.stats_mode = 2;

do_both_animals=0;
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
        perm_mode=41.6013111;       
        script_name = 'plots_CFC';
        freqband_stats = freqband_stats_alpha;
        force_freqband_stats = 1;
end

% [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    
%exclude_clipping = 1;  warning('Including clip!');
exclude_60 = 0;

permdat2pls = 1;            % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    do_diff = 1;            % Actually, don't normalize by null distribution
    sort_pls = 0;
perm2pls = 0;               % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
    perm2pls_dophi = 0;
    
plot_SFC_spectra=0;


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

%%
% inds = [1:length(groupls)];
% inds = [1,2,3,4,9,10];
% names = {'SchA Pref','SchA NonPref','SchB Pref','SchB NonPref','Pref Sch','NonPref Sch'};

% inds = [1,2,5];
% names = {'SchA Pref - NonPref','SchB Pref - NonPref','Pref Sch - NonPref Sch'};
% names = {'Ctg1 vs 2','Ctg3 vs 4','SchA vs B'};
% 
inds = [1:3];
names = {'Ctg 1 vs 2','Ctg 3 vs 4','Sch A vs B'};
% names = {'SchA Pref - NonPref','SchB Pref - NonPref','Pref Sch - NonPref Sch'};

[groupls(inds).legend] = names{:};
[groupld(inds).legend] = names{:};
[groupos(inds).legend] = names{:};
[groupod(inds).legend] = names{:};


i=0;
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupls(inds) ],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupld(inds) ],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupos(inds) ],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupod(inds) ],[],opts_PM3Dcs,[]);

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
opts_PPS_temp.do_bh0=1;
opts_PPS_temp.do_phi=perm2pls_dophi;

plot_scatter_both_animals = 1;
include_SchAvsB = 1;

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
    if include_SchAvsB;
        if i ~=3
            sp_ls = [sp_ls PPS_temp.q4_perm(:)];
        else
            sp_ls = [sp_ls PPS_temp.q4_permp(:) PPS_temp.q4_permm(:)];
        end
    end
    clear PPS_temp
end

N_criteria = 4;
gr_template = Grp;
gr_template.criteria = [2*ones(1,N_criteria)]; gr_template.criteria_alt = []; gr_template.criteria_sfc = []; gr_template.ctgs = 1;
gr_template.xlims_desired = [0 120];


clear group
i=0;
% Ctg A Sensitive - Pref vs non-preferred
% i=i+1; group(i) = gr_template; group(i).criteria=[1 0 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch A';
% i=i+1; group(i) = gr_template; group(i).criteria=[0 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch B';
% i=i+1; group(i) = gr_template; group(i).criteria=[1 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Both';
i=i+1; group(i) = gr_template; group(i).criteria=[2 2 1 0]; group(i).ctgs = [1]; group(i).legend = 'Sch+';
i=i+1; group(i) = gr_template; group(i).criteria=[2 2 0 1]; group(i).ctgs = [1]; group(i).legend = 'Sch-';

for i = 1:length(group); [group(i).cells, group(i).numcells]= get_grouped_cells(group(i),[sp_ls]); end

if plot_scatter_both_animals
    figure; plot(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),'kx','LineWidth',2); hold on; plot_scattergroups(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),group,~bad_curr); xlabel('FFC Ctg1-Ctg2'); ylabel('FFC Ctg3-Ctg4')
else
    figure; plot(pls_curr(:,1),pls_curr(:,2),'k.'); hold on; plot_scattergroups(pls_curr(:,1),pls_curr(:,2),group); xlabel('FFC Ctg1-Ctg2'); ylabel('FFC Ctg3-Ctg4')
end
title('Monkey L Sample');

%% Monkey L Delay
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
    if include_SchAvsB;
        if i ~=3
            sp_ls = [sp_ls PPS_temp.q4_perm(:)];
        else
            sp_ls = [sp_ls PPS_temp.q4_permp(:) PPS_temp.q4_permm(:)];
        end
    end
    clear PPS_temp
end

N_criteria = 4;
gr_template = Grp;
gr_template.criteria = [2*ones(1,N_criteria)]; gr_template.criteria_alt = []; gr_template.criteria_sfc = []; gr_template.ctgs = 1;
gr_template.xlims_desired = [0 120];


clear group
i=0;
% Ctg A Sensitive - Pref vs non-preferred
% i=i+1; group(i) = gr_template; group(i).criteria=[1 0 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch A';
% i=i+1; group(i) = gr_template; group(i).criteria=[0 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch B';
% i=i+1; group(i) = gr_template; group(i).criteria=[1 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Both';
i=i+1; group(i) = gr_template; group(i).criteria=[2 2 1 0]; group(i).ctgs = [1]; group(i).legend = 'Sch+';
i=i+1; group(i) = gr_template; group(i).criteria=[2 2 0 1]; group(i).ctgs = [1]; group(i).legend = 'Sch-';

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
    if include_SchAvsB;
        if i ~=3
            sp_ls = [sp_ls PPS_temp.q4_perm(:)];
        else
            sp_ls = [sp_ls PPS_temp.q4_permp(:) PPS_temp.q4_permm(:)];
        end
    end
    clear PPS_temp
end

N_criteria = 4;
gr_template = Grp;
gr_template.criteria = [2*ones(1,N_criteria)]; gr_template.criteria_alt = []; gr_template.criteria_sfc = []; gr_template.ctgs = 1;
gr_template.xlims_desired = [0 120];


clear group
i=0;
% Ctg A Sensitive - Pref vs non-preferred
% i=i+1; group(i) = gr_template; group(i).criteria=[1 0 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch A';
% i=i+1; group(i) = gr_template; group(i).criteria=[0 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch B';
% i=i+1; group(i) = gr_template; group(i).criteria=[1 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Both';
i=i+1; group(i) = gr_template; group(i).criteria=[2 2 1 0]; group(i).ctgs = [1]; group(i).legend = 'Sch+';
i=i+1; group(i) = gr_template; group(i).criteria=[2 2 0 1]; group(i).ctgs = [1]; group(i).legend = 'Sch-';

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
    if include_SchAvsB;
        if i ~=3
            sp_ls = [sp_ls PPS_temp.q4_perm(:)];
        else
            sp_ls = [sp_ls PPS_temp.q4_permp(:) PPS_temp.q4_permm(:)];
        end
    end
    clear PPS_temp
end

N_criteria = 4;
gr_template = Grp;
gr_template.criteria = [2*ones(1,N_criteria)]; gr_template.criteria_alt = []; gr_template.criteria_sfc = []; gr_template.ctgs = 1;
gr_template.xlims_desired = [0 120];


clear group
i=0;
% Ctg A Sensitive - Pref vs non-preferred
% i=i+1; group(i) = gr_template; group(i).criteria=[1 0 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch A';
% i=i+1; group(i) = gr_template; group(i).criteria=[0 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Sch B';
% i=i+1; group(i) = gr_template; group(i).criteria=[1 1 2*ones(1,N_criteria-2)]; group(i).ctgs = [1]; group(i).legend = 'Both';
i=i+1; group(i) = gr_template; group(i).criteria=[2 2 1 0]; group(i).ctgs = [1]; group(i).legend = 'Sch+';
i=i+1; group(i) = gr_template; group(i).criteria=[2 2 0 1]; group(i).ctgs = [1]; group(i).legend = 'Sch-';

for i = 1:length(group); [group(i).cells, group(i).numcells]= get_grouped_cells(group(i),[sp_ls]); end

if plot_scatter_both_animals
    figure; plot(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),'kx','LineWidth',2); hold on; plot_scattergroups(pls_curr(~bad_curr,1),pls_curr(~bad_curr,2),group,~bad_curr); xlabel('FFC Ctg1-Ctg2'); ylabel('FFC Ctg3-Ctg4')
else
    figure; plot(pls_curr(:,1),pls_curr(:,2),'k.'); hold on; plot_scattergroups(pls_curr(:,1),pls_curr(:,2),group); xlabel('FFC Ctg1-Ctg2'); ylabel('FFC Ctg3-Ctg4')
end
title('Monkey O Delay');

% Fig 7 - Ctgs Scatterplots
