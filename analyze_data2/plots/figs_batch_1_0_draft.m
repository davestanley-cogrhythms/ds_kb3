
% % % % Random leftover things
% 
% % Old save command
% saveas(gcf,savenames{i},savetype);
% 
% % For loop for saving figs
if ~exist('currfname'); currfname = 'figs_batch_2_0_ctgs'; end
if ~exist('currfigname'); currfigname = 'Fig6'; end
savenames={'fig1','fig2','fig3','fig4','fig5','fig6','fig7','fig8','fig9','fig10','fig11','fig12','fig13','fig14','fig15','fig16','fig17','fig18'};
mydate = datestr(datenum(date),'yy/mm/dd'); mydate = strrep(mydate,'/','');
c=clock;
sp = ['d' mydate '_t' num2str(c(4),'%10.2d') '' num2str(c(5),'%10.2d') '' num2str(round(c(6)),'%10.2d')];
sp = [sp '__' currfname '_' currfigname];
mkdir(fullfile('.',sp));
for i=[1:gcf]
    figure(i); %ylim([0 15])
    set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',fullfile('.',sp,savenames{i}))
end
% 
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
    
%% Fig 0 - Supp - FFC Alpha vs Beta peaks regression

    clearvars -except wrkspc_buffer fv
    vars_pull(fv);
    currfigname = 'Fig0';

    sfc_mode=22.4013100;
    perm_mode=22.4013110;
        [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);

    exclude_60 = 0;
    plot_regress_alpha_vs_beta = 1;

    
    % % Stage & Animals
    curr_stage_sfc=3;
        curr_stage_badclipping = curr_stage_sfc;
    excludeL = 0; excludeO = 1; % Animals
        
    % Monkey L
    excludeL = 0; excludeO = 1; 
    push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
    
    % Monkey O
    excludeL = 1; excludeO = 0; 
    push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
    
%% Fig 0b - Supp - SFC beta peak vs FFC peaks regression

    clearvars -except wrkspc_buffer fv
    vars_pull(fv);
    currfigname = 'Fig0b';

    sfc_mode=22.4013100;
    perm_mode=22.4013110;
        [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);

    exclude_60 = 0;
    exclude_lowfiring = 1;
    
    
    regress_sfc_vs_ffc2 = 1;

    
    % % Stage & Animals
    curr_stage_sfc=3;
        curr_stage_badclipping = curr_stage_sfc;
    excludeL = 0; excludeO = 1; % Animals
        
    % Monkey L
    excludeL = 0; excludeO = 1; 
    push_workspace; out = script2func('plots_preferred_unitpairs'); wrkspc_buffer = out.wrkspc_buffer; 
    
    % Monkey O
    excludeL = 1; excludeO = 0; 
    push_workspace; out = script2func('plots_preferred_unitpairs'); wrkspc_buffer = out.wrkspc_buffer; 
    
    

%% Fig 1A - Sample vs Delay - FFC - Spectrum
clearvars -except wrkspc_buffer fv
vars_pull(fv); close all
currfigname = 'Fig1A';

sfc_mode=22.4013100;
perm_mode=22.4013110;
boot_mode = 22.401312;  % Not used; this is all we have right now
    [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    
    
N_criteria = 3;
gr2 = Grp;
gr2.criteria = [2*ones(1,N_criteria)]; gr2.criteria_alt = [2*ones(1,N_criteria)]; gr2.ctgs = 11;
group=gr2;

% % Stage
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;

excludeL = 0; excludeO = 1; % Animals
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
groupls=out.group; clear out
groupls.legend='Sample';

excludeL = 1; excludeO = 0; % Animals
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
groupos=out.group; clear out
groupos.legend='Sample';

% % Stage
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;

excludeL = 0; excludeO = 1; % Animals
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
groupld=out.group; clear out
groupld.legend='Delay';

excludeL = 1; excludeO = 0; % Animals
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
groupod=out.group; clear out
groupod.legend='Delay';

i=0;
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupls groupld],[],opts_PM3Dcs,['bgbg']);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupos groupod],[],opts_PM3Dcs,['bgbg']);




%% Fig 1A1 - Sample vs Delay - FFC - Phi
clearvars -except wrkspc_buffer fv
vars_pull(fv); close all
currfigname = 'Fig1A1';

sfc_mode=22.4013100;
perm_mode=22.4013110;
boot_mode = 22.401312;  % Not used; this is all we have right now
    [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    
plotmode = 8;               % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    
N_criteria = 3;
gr2 = Grp;
gr2.criteria = [2*ones(1,N_criteria)]; gr2.criteria_alt = [2*ones(1,N_criteria)]; gr2.ctgs = 11;
group=gr2;


freqband_stats = freqband_stats_alpha;

% % Stage
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;

excludeL = 0; excludeO = 1; % Animals
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
groupls=out.group; clear out
groupls.legend='Sample';

excludeL = 1; excludeO = 0; % Animals
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
groupos=out.group; clear out
groupos.legend='Sample';

% % Stage
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;

excludeL = 0; excludeO = 1; % Animals
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
groupld=out.group; clear out
groupld.legend='Delay';

excludeL = 1; excludeO = 0; % Animals
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
groupod=out.group; clear out
groupod.legend='Delay';

grouplsA=groupls;
groupldA=groupld;
grouposA=groupos;
groupodA=groupod;


freqband_stats = freqband_stats_lbeta;

% % Stage
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;

excludeL = 0; excludeO = 1; % Animals
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
groupls=out.group; clear out
groupls.legend='Sample';

excludeL = 1; excludeO = 0; % Animals
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
groupos=out.group; clear out
groupos.legend='Sample';

% % Stage
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;

excludeL = 0; excludeO = 1; % Animals
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
groupld=out.group; clear out
groupld.legend='Delay';

excludeL = 1; excludeO = 0; % Animals
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
groupod=out.group; clear out
groupod.legend='Delay';

grouplsB=groupls;
groupldB=groupld;
grouposB=groupos;
groupodB=groupod;

x = groupldA.datastats; y = groupldB.datastats;
% figure; plot(abs(x),angle(x),'b.'); hold on; plot(abs(y),angle(y),'r.');
% figure; [p,n] = ecdf(angle(x)); plot(n,p,'b'); hold on; [p,n] = ecdf(angle(y)); plot(n,p,'r');
figure; hold on; plot(x,'.'); plot(y,'.'); plot(y,'r.'); 
    plot(mean(x),'x','MarkerSize',50,'LineWidth',10); viscircles([0,0],[abs(mean(x))],'EdgeColor','b');
    plot(mean(y),'rx','MarkerSize',50,'LineWidth',10); viscircles([0,0],[abs(mean(y))],'EdgeColor','r');

x = groupodA.datastats; y = groupodB.datastats;
% figure; plot(abs(x),angle(x),'b.'); hold on; plot(abs(y),angle(y),'r.');
% figure; [p,n] = ecdf(angle(x)); plot(n,p,'b'); hold on; [p,n] = ecdf(angle(y)); plot(n,p,'r');
figure; hold on; plot(x,'.'); plot(y,'.'); plot(y,'r.'); 
    plot(mean(x),'x','MarkerSize',50,'LineWidth',10); viscircles([0,0],[abs(mean(x))],'EdgeColor','b');
    plot(mean(y),'rx','MarkerSize',50,'LineWidth',10); viscircles([0,0],[abs(mean(y))],'EdgeColor','r');

opts_PSC.hmask = blkdiag(ones(2),ones(2),ones(2));
i=i+1; [h1,hstats,pstats,out] = plot_phistats_custstruct([ groupldA(1) groupodA(1) ...
    grouplsB(1) grouposB(1)  groupldB(1) groupodB(1)  ...
    ],opts_PSC)



plot_SFC_spectra = 1;


% % Stage
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;

plotmode = 1;
excludeL = 0; excludeO = 1; % Animals
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
groupld=out.group; clear out
groupld.legend='E(|FFC|)';

plotmode = 6;
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
groupld=out.group; clear out
groupld.legend='|E(FFC)|';

excludeL = 1; excludeO = 0; % Animals
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
groupod=out.group; clear out
groupod.legend='Delay';

%% Fig 1B - Sample vs Delay All - SFC - Spectrum
clearvars -except wrkspc_buffer fv
vars_pull(fv); close all
currfigname = 'Fig1B';

sfc_mode=2.2015100;
perm_mode = 2.201511;
boot_mode = 2.201512;
    [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    

% % Group Setup
symmetry_mode=3;
[gr,mask] = group_setup_FR2D;       % Generate grouping structure
[gr.all_cells, gr.all_sens, gr.all_sens4, gr.single, gr.single_sens] = group_setup_other;

group = gr.all_cells;
group.ctgs = 11;
    
% % Stage
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;

excludeL = 0; excludeO = 1; % Animals
push_workspace; out = script2func('plots_preferred'); wrkspc_buffer = out.wrkspc_buffer; 
groupls=out.group; clear out
groupls.legend='Sample';

excludeL = 1; excludeO = 0; % Animals
push_workspace; out = script2func('plots_preferred'); wrkspc_buffer = out.wrkspc_buffer; 
groupos=out.group; clear out
groupos.legend='Sample';

% % Stage
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;

excludeL = 0; excludeO = 1; % Animals
push_workspace; out = script2func('plots_preferred'); wrkspc_buffer = out.wrkspc_buffer; 
groupld=out.group; clear out
groupld.legend='Delay';

excludeL = 1; excludeO = 0; % Animals
push_workspace; out = script2func('plots_preferred'); wrkspc_buffer = out.wrkspc_buffer; 
groupod=out.group; clear out
groupod.legend='Delay';

i=0;
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupls groupld],[],opts_PM3Dcs,['bgbg']); set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupos groupod],[],opts_PM3Dcs,['bgbg']); set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})




%% Fig 1C - Sample vs Delay All - PSD - Spectrum
clearvars -except wrkspc_buffer fv
vars_pull(fv); close all
currfigname = 'Fig1C';

psd_mode=41.60131;

[~, mode_subgroups] = decode_sfc_mode(psd_mode);
[fname_suffix] = build_sfcmode(psd_mode, mode_subgroups);


% % Group setup
group_mode = 1;
N_criteria = 3;
gr2 = Grp;
gr2.criteria = [2*ones(1,N_criteria)]; gr2.criteria_alt = [2*ones(1,N_criteria)]; gr2.ctgs = 11;
group=gr2;

% % Stage
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;

excludeL = 0; excludeO = 1; % Animals
push_workspace; out = script2func('plots_CFC'); wrkspc_buffer = out.wrkspc_buffer; 
groupls=out.group; clear out
groupls.legend='Sample';

excludeL = 1; excludeO = 0; % Animals
push_workspace; out = script2func('plots_CFC'); wrkspc_buffer = out.wrkspc_buffer; 
groupos=out.group; clear out
groupos.legend='Sample';

% % Stage
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;

excludeL = 0; excludeO = 1; % Animals
push_workspace; out = script2func('plots_CFC'); wrkspc_buffer = out.wrkspc_buffer; 
groupld=out.group; clear out
groupld.legend='Delay';

excludeL = 1; excludeO = 0; % Animals
push_workspace; out = script2func('plots_CFC'); wrkspc_buffer = out.wrkspc_buffer; 
groupod=out.group; clear out
groupod.legend='Delay';

i=0;
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupls groupld],[],opts_PM3Dcs,['bgbg']); set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupos groupod],[],opts_PM3Dcs,['bgbg']); set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})

%% Fig 1D - Sample vs Delay - Spectrograms (SFC, LFP, spiking)

clearvars -except wrkspc_buffer fv
vars_pull(fv);
currfigname = 'Fig1D';

plot_spectrogram = 1;

% Decide what spectrograms to plot (spikes/sfc/lfp)
myspecmode = 3; %1-SFC; 2-LFP; 3-spikes
switch myspecmode
    case 1
        spect_mode = 3.20101;
        plotmode=1;
    case 2
        spect_mode = 8.00001;
        plotmode=1;
        exclude_lowfiring = 0;
    case 3
        spect_mode = 3.20101;
        plotmode=5;
        exclude_clipping = 0;
        exclude_60 = 0;
        exclude_lowfiring = 0;
end


% % Group Setup
symmetry_mode=3;
[gr,mask] = group_setup_FR2D;       % Generate grouping structure
[gr.all_cells, gr.all_sens, gr.all_sens4, gr.single, gr.single_sens] = group_setup_other;

group = gr.all_cells;
group.ctgs = 11;

excludeL = 0; excludeO = 1; % Animals
push_workspace; out = script2func('plots_preferred'); wrkspc_buffer = out.wrkspc_buffer; 
group1a=out.group; clear out
group1a.legend='Monkey L';

excludeL = 1; excludeO = 0; % Animals
push_workspace; out = script2func('plots_preferred'); wrkspc_buffer = out.wrkspc_buffer; 
group1b=out.group; clear out
group1b.legend='Monkey O';



%% Fig 2A - Sch AvsB
clearvars -except wrkspc_buffer fv
vars_pull(fv); %close all
currfigname = 'Fig2A';

sfc_mode=22.4013100;
perm_mode=22.4013110;
boot_mode = 22.401312;  % Not used; this is all we have right now
    [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    

%freqband_stats = fv.freqband_stats_alpha;
freqband_stats = fv.freqband_stats_lbeta;

% Setup gr2
N_criteria = 3;
gr2 = Grp;
gr2.criteria = [2*ones(1,N_criteria)]; gr2.criteria_alt = [2*ones(1,N_criteria)]; gr2.ctgs = 1;

clear group
i=0;
i=i+1; group(i) = gr2; group(i).ctgs = [9]; group(i).legend = 'SchA';
i=i+1; group(i) = gr2; group(i).ctgs = [10]; group(i).legend = 'SchB';     

opts_PSC.hmask = blkdiag(ones(2),ones(2),ones(2),ones(2));

% % Monkey L
excludeL = 0; 
excludeO = 1; 

curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func('plots_ffc_vs_Sch');
out.group(1).legend = 'Sample SchA';
out.group(2).legend = 'Sample SchB';
group1=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out

curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func('plots_ffc_vs_Sch');
out.group(1).legend = 'Delay SchA';
out.group(2).legend = 'Delay SchB';
group2=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out

% % Monkey O
excludeL = 1; 
excludeO = 0; 

curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func('plots_ffc_vs_Sch');
out.group(1).legend = 'Sample SchA';
out.group(2).legend = 'Sample SchB';
group3=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out

curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func('plots_ffc_vs_Sch');
out.group(1).legend = 'Delay SchA';
out.group(2).legend = 'Delay SchB';
group4=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out


% Plot figures
i=0;
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[group1 group2],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[group3 group4],[],opts_PM3Dcs,[]);
i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group1 group2 group3 group4],opts_PSC);

% Save figures
for i=1:3
%     figure(i);
%     set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})
end

%% Fig 2B - Sch AvB - FFC - Pref vs non-preferrred Spectrum
clearvars -except wrkspc_buffer fv
vars_pull(fv); close all
currfigname = 'Fig2B';

sfc_mode=22.4013100;
perm_mode=22.4013110;
boot_mode = 22.401312;  % Not used; this is all we have right now
    [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    

% % % Setup group
[gr,mask] = group_setup_FR2D;       % Generate grouping structure

clear group1; N1 = length(gr.all); group1(1:N1) = Grp_SD; for i = 1:N1; group1(i) = group1(i).inheritObj(gr.all(i)); end
clear group2; N2 = length(gr.all); group2(1:N2) = Grp_SD; for i = 1:N2; group2(i) = group2(i).inheritObj(gr.all(i)); end

group2 = group_swap.alt(group2);
groupLR = group_OUTER(group1,group2);

% (Groupmode = 11)
clear group                             
i=0;
i=i+1;group(i) = group_merge(groupLR([11,12,15,3,9,6,8,14,2,5])); group(i).ctgs = [9;9;9;9;9;10;10;10;10;10]; group(i).legend=('OR Preferred Sch');
i=i+1;group(i) = group_merge(groupLR([11,12,15,3,9,6,8,14,2,5])); group(i).ctgs = [10;10;10;10;10;9;9;9;9;9]; group(i).legend=('OR Nonpreferred Sch');

opts_PSC.hmask = blkdiag(ones(2),ones(2),ones(2),ones(2));

% % Stage
curr_stage=2;
curr_stage_alt=3;
curr_stage_sfc=2;
curr_stage_badclipping = curr_stage_sfc;

% Animals
excludeL = 0; 
excludeO = 1; 

push_workspace
out = script2func('plots_preferred_unitpairs');
group1a=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out
group1a(1).legend = 'Sample Preferred';
group1a(2).legend = 'Sample Non-Preferred';

% Animals
excludeL = 1; 
excludeO = 0; 

push_workspace
out = script2func('plots_preferred_unitpairs');
group1b=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out
group1b(1).legend = 'Sample Preferred';
group1b(2).legend = 'Sample Non-Preferred';


% % Stage
curr_stage=3;
curr_stage_alt=3;
curr_stage_sfc=3;
curr_stage_badclipping = curr_stage_sfc;

% Animals
excludeL = 0; 
excludeO = 1; 

push_workspace
out = script2func('plots_preferred_unitpairs');
group2a=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out
group2a(1).legend = 'Delay Preferred';
group2a(2).legend = 'Delay Non-Preferred';

% Animals
excludeL = 1; 
excludeO = 0; 

push_workspace
out = script2func('plots_preferred_unitpairs');
group2b=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out
group2b(1).legend = 'Delay Preferred';
group2b(2).legend = 'Delay Non-Preferred';

% Plot figures
i=0;
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[group1a group2a],[],opts_PM3Dcs,[]); set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[group1b group2b],[],opts_PM3Dcs,[]); set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})
i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group1a group2a group1b group2b],opts_PSC);

% Save figures
for i=1:3
    figure(i);
    set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})
end


%figure; [h1, out.PM3Dcs] = plot_matrix3D_custstruct([],[group1a group2a],[],[],['bgbg']);
% figure; [h1, out.PM3Dcs] = plot_matrix3D_custstruct([],[group1a group2a],[],[],[0 1 0; 0 0.8 0; 0 0 1; 0 0 0.8]);
%figure; [h1, out.PM3Dcs] = plot_matrix3D_custstruct([],[group1b group2b]);

%% Fig 2C - Sch AvB - FFC - Scatterplot

clearvars -except wrkspc_buffer fv
vars_pull(fv); close all
currfigname = 'Fig2C';

sfc_mode=22.4013100;
perm_mode=22.4013110;
boot_mode = 22.401312;  % Not used; this is all we have right now
    [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    


exclude_60 = 0;
excludeL = 0; 
excludeO = 0; 

sens_chosen_condition = 1;
sens_ctg1_name='Sch A';
sens_ctg2_name='Sch B';

plot_Sch_sens = 1;

% % Stage
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;

freqband_stats = fv.freqband_stats_alpha;
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
perm1a=out.PPS.q4_perm; bad_any = out.PPS.bad_any; if sum(bad_any) > 0; warning([num2str(bad_any) ' bad cells detected']); end;
clear out

freqband_stats = fv.freqband_stats_lbeta;
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
perm1b=out.PPS.q4_perm; bad_any = out.PPS.bad_any; if sum(bad_any) > 0; warning([num2str(bad_any) ' bad cells detected']); end;
clear out


% % Stage
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;

freqband_stats = fv.freqband_stats_alpha;
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
perm2a=out.PPS.q4_perm; bad_any = out.PPS.bad_any; if sum(bad_any) > 0; warning([num2str(bad_any) ' bad cells detected']); end;
clear out

freqband_stats = fv.freqband_stats_lbeta;
push_workspace; out = script2func('plots_ffc_vs_Sch'); wrkspc_buffer = out.wrkspc_buffer; 
perm2b=out.PPS.q4_perm; bad_any = out.PPS.bad_any; if sum(bad_any) > 0; warning([num2str(bad_any) ' bad cells detected']); end;
clear out

% Save figures
for i=1:4
%     figure(i);
%     set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})
end

% Make pie graph
% figure; plot_venn(perm1b, perm2b,perm1a,{'A','B','C'},0,1)
% Make bar graph
i=i+1;
data = [perm1a(:), perm2a(:), perm1b(:), perm2b(:)];

warning('this is not quite correct b/c includes bad files');
databs = bootstrp(150,@sum,data);
%databs=databs/length(data)*100;     % Change to percent
fign; bh=boxplot(databs(:,:),'boxstyle','outline','colors','k','symbol','r+');
set(bh(:),'linewidth',4);
set(gca,'XTickLabel',{''},'Box','off');
ylabel('N electrode pairs'); 
set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})

% labels = repmat({'Alpha S','Alpha D','Beta S','Beta D'},150,1);
% databs=databs(:);
% labels=labels(:);
% fign; boxplot(databs(:),labels(:))


%% Fig 2C2 - Sch AvB - All - PermuteSpect

clearvars -except wrkspc_buffer fv
vars_pull(fv);
currfigname = 'Fig2C2';



ffc_sfc_psd_mode = 1;
switch ffc_sfc_psd_mode
    case 1
        sfc_mode=22.4013100;    % AvB
        perm_mode=22.4013110;
        %sfc_mode=22.4313101;    % RT
        %perm_mode=22.4313111;
        script_name = 'plots_ffc_vs_Sch';
    case 2
        sfc_mode=2.2015100;
        perm_mode=2.2015110;
        script_name = 'plots_preferred';
    case 3
        sfc_mode=41.6013100;
        perm_mode=41.6013110;       % Incomplete
        script_name = 'plots_CFC';
end

[~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    

exclude_60 = 0;

perm2pls = 1;               % Instead of showing raw SFC values, show % cells successfully passing permutation test
        perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
        perm2pls_dophi = 0;

plot_SFC_spectra=0;

% Choose animals & stage
excludeL = 0; 
excludeO = 1; 
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer; out.group.legend = 'L Sample';
groupls = out.group;
clear out

% Choose animals & stage
excludeL = 0; 
excludeO = 1; 
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer; out.group.legend = 'L Delay';
groupld = out.group;
clear out

% Choose animals & stage
excludeL = 1; 
excludeO = 0; 
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer; out.group.legend = 'O Sample';
groupos = out.group;
clear out

% Choose animals & stage
excludeL = 1; 
excludeO = 0; 
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer; out.group.legend = 'O Delay';
groupod = out.group;
clear out

opts_PM3Dcs.do_sgolay=1;
opts_PM3Dcs.stats_mode = 2;
i=0;
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupls groupld groupos groupod],[],opts_PM3Dcs,[]); xlim([0 200])


% Save figures
% for i=1:1
%     figure(i);
%     set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})
% end



%% Fig2D - Sch AvB - SFC - Pref vs non-preferrred Spect OR AvB
clearvars -except wrkspc_buffer fv
vars_pull(fv);
currfigname = 'Fig2D';

do_AvB = 0;

sfc_mode=2.2015100;
perm_mode = 2.201511;
boot_mode = 2.201512;
    [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);


% % Group Setup
symmetry_mode=3;
[gr,mask] = group_setup_FR2D;       % Generate grouping structure
[gr.all_cells, gr.all_sens, gr.all_sens4, gr.single, gr.single_sens] = group_setup_other;

% Scored by Right minus wRong
%gr.all = normalize_diffs(gr.all,'diff');
gr.rmr = score2grouping(gr.all,[gr.all.diff]);                     % Don't use Mask for rmr
% gr.rmr = score2grouping(gr.all(mask),[gr.all(mask).diff]);         % Use Mask for rmr
% gr_netpref = score2grouping(gr.all,[gr.all.netprefA]); gr_netpref.ctgs = []

gr = group_enumerate_sens(gr,mask,symmetry_mode);     % Enumerate grouping structure across various "sensitivity tests"
group = gr.sens_sch;
if do_AvB; group = gr.AvB_everything;end
    
% % Stage
curr_stage=2;
curr_stage_alt=3;
curr_stage_sfc=2;
curr_stage_badclipping = curr_stage_sfc;

% Animals
excludeL = 0; 
excludeO = 1; 

push_workspace
out = script2func('plots_preferred');
group1a=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out
group1a(1).legend = 'Sample Preferred';
group1a(2).legend = 'Sample Non-Preferred';
if do_AvB
    group1a(1).legend = 'Sample Sch A';
    group1a(2).legend = 'Sample Sch B'; end


% Animals
excludeL = 1; 
excludeO = 0; 

push_workspace
out = script2func('plots_preferred');
group1b=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out
group1b(1).legend = 'Sample Preferred';
group1b(2).legend = 'Sample Non-Preferred';
if do_AvB
    group1b(1).legend = 'Delay Sch A';
    group1b(2).legend = 'Delay Sch B'; end

% % Stage
curr_stage=3;
curr_stage_alt=3;
curr_stage_sfc=3;
curr_stage_badclipping = curr_stage_sfc;

% Animals
excludeL = 0; 
excludeO = 1; 

push_workspace
out = script2func('plots_preferred');
group2a=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out
group2a(1).legend = 'Delay Preferred';
group2a(2).legend = 'Delay Non-Preferred';
if do_AvB
    group2a(1).legend = 'Sample Sch A';
    group2a(2).legend = 'Sample Sch B'; end

% Animals
excludeL = 1; 
excludeO = 0; 

push_workspace
out = script2func('plots_preferred');
group2b=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out
group2b(1).legend = 'Delay Preferred';
group2b(2).legend = 'Delay Non-Preferred';
if do_AvB
    group2b(1).legend = 'Delay Sch A';
    group2b(2).legend = 'Delay Sch B'; end

i=0;
opts_PSC.hmask = blkdiag(true(2),true(2),true(2),true(2));
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[group1a group2a],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[group1b group2b],[],opts_PM3Dcs,[]);
i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group1a group2a group1b group2b],opts_PSC);



%% Fig 3A1 - RT - Spectra
clearvars -except wrkspc_buffer fv
vars_pull(fv); %close all
currfigname = 'Fig3A1';

gr_submode = 1;
opts_PSC.hmask = blkdiag(ones(2),ones(2),ones(2),ones(2));

ffc_sfc_psd_mode = 1;
switch ffc_sfc_psd_mode
    case 1
        sfc_mode=22.4213101;    % RT
        %perm_mode=22.4213111;  % Not yet implemented
        script_name = 'plots_ffc_vs_Sch';
        N_criteria=3;
        %freqband_stats = fv.freqband_stats_lbeta;
        freqband_stats = fv.freqband_stats_alpha;
    case 2
        sfc_mode=2.4215100;
        perm_mode=2.4215110;
        script_name = 'plots_preferred';
        N_criteria=5;   % For Roy
        freqband_stats = fv.freqband_stats_hbeta;
    case 3
        sfc_mode=41.6213101;
        perm_mode=41.6213111;
        script_name = 'plots_CFC';
        N_criteria=3;
        freqband_stats = fv.freqband_stats_alpha;
        %freqband_stats = fv.freqband_stats_lbeta;
end

[~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);

    
group_override = get_group_overrides(fname_suffix_sfc, N_criteria, gr_submode);
if ~isempty(group_override); group = group_override; end


% % Monkey L
excludeL = 0; 
excludeO = 1; 

exclude_60 = 0;

curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func(script_name);
groupls=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out

curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func(script_name);
groupld=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out

% % Monkey O
excludeL = 1; 
excludeO = 0; 

curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func(script_name);
groupos=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out

curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func(script_name);
groupod=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out


append_samp = @(s) ['Samp. ',s];
append_delay = @(s) ['Delay ',s];
groupls=myfunc(groupls,'legend',append_samp)
groupld=myfunc(groupld,'legend',append_delay)
groupos=myfunc(groupos,'legend',append_samp)
groupod=myfunc(groupod,'legend',append_delay)


% Plot figures
i=0;
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupls groupld],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupos groupod],[],opts_PM3Dcs,[]);
i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([groupls groupld groupos groupod],opts_PSC);

% Save figures
for i=1:3
%     figure(i);
%     set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})
end

%% Fig 3A2 - RT - All - PermuteSpect

clearvars -except wrkspc_buffer fv
vars_pull(fv);
currfigname = 'Fig3A2';

opts_PM3Dcs.do_sgolay = 1;

ffc_sfc_psd_mode = 3;
switch ffc_sfc_psd_mode
    case 1
        sfc_mode=22.4213101;    % RT
        perm_mode=22.4213111;
        script_name = 'plots_ffc_vs_Sch';
    case 2
        sfc_mode=2.4215100;
        perm_mode=2.4215110;
        script_name = 'plots_preferred';
    case 3
        sfc_mode=41.6213101;
        perm_mode=41.6213111;     
        script_name = 'plots_CFC';
end

[~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    
exclude_60 = 0;

perm2pls = 1;               % Instead of showing raw SFC values, show % cells successfully passing permutation test
        perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
        perm2pls_dophi = 0;

plot_SFC_spectra=0;

% Choose animals & stage
excludeL = 0; 
excludeO = 1; 
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
groupls = out.group;
clear out

% Choose animals & stage
excludeL = 0; 
excludeO = 1; 
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
groupld = out.group;
clear out

% Choose animals & stage
excludeL = 1; 
excludeO = 0; 
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
groupos = out.group;
clear out

% Choose animals & stage
excludeL = 1; 
excludeO = 0; 
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
groupod = out.group;
clear out

groupls(1).legend = 'Sample Congruent RT';
groupld(1).legend = 'Delay Congruent RT';
groupos(1).legend = 'Sample Congruent RT';
groupod(1).legend = 'Delay Congruent RT';

i=0;
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupls groupld],[],opts_PM3Dcs,[]);
% i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupos groupod],[],opts_PM3Dcs,[]);
% i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[],[],opts_PM3Dcs,[]);


% Save figures
% for i=1:2; figure(i); set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i}); end


%% Fig 3B1 - RT2 - Spectra


clearvars -except wrkspc_buffer fv
vars_pull(fv); %close all
currfigname = 'Fig3B1';

gr_submode = 1:2;

ffc_sfc_psd_mode = 1;
switch ffc_sfc_psd_mode
    case 1
        sfc_mode=22.4313101;    % RT
        perm_mode=22.4313111;
        script_name = 'plots_ffc_vs_Sch';
        N_criteria=3;
        freqband_stats = fv.freqband_stats_lbeta;
    case 2
        sfc_mode=2.4315100;
        perm_mode=2.4315110;
        script_name = 'plots_preferred';
        N_criteria=5;   % For Roy
        freqband_stats = fv.freqband_stats_hbeta;
    case 3
        sfc_mode=41.6313101;
        perm_mode=41.6313111;
        script_name = 'plots_CFC';
        N_criteria=3;
        freqband_stats = fv.freqband_stats_alpha;
        %freqband_stats = fv.freqband_stats_lbeta;
end

[~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);

    
group_override = get_group_overrides(fname_suffix_sfc, N_criteria, gr_submode);
if ~isempty(group_override); group = group_override; end


%opts_PSC.hmask = blkdiag(ones(2),ones(2),ones(2),ones(2));

% % Monkey L
excludeL = 0; 
excludeO = 1; 

exclude_60 = 0;

curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func(script_name);
groupls=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out

curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func(script_name);
groupld=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out

% % Monkey O
excludeL = 1; 
excludeO = 0; 

curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func(script_name);
groupos=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out

curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func(script_name);
groupod=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out

cind = 1:2;
iind = 3:4;

switch ffc_sfc_psd_mode;
    case 1; units= 'FFC'; scaling=-1;
    case 3; units= 'PSD'; scaling=1;
end

%%
i=0;
% % Plot figures SAMPLE
opts_PSC.hmask = blkdiag(ones(2),ones(2),ones(2),ones(2));
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupls(cind)],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupls(iind)],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupos(cind)],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupos(iind)],[],opts_PM3Dcs,[]);
i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([groupls groupos],opts_PSC);
% 
% %%

% % Plot Diffs SAMPLE
opts_PSC.hmask = blkdiag(ones(2),ones(2));
clear group; group(1:4) = Grp; group(1).freqband_stats = freqband_stats;
group(1).datastats = scaling*(groupls(1).datastats - groupls(2).datastats)./mean(groupls(2).datastats)*100; group(1).legend=['Congruent'];
group(2).datastats = scaling*(groupls(3).datastats - groupls(4).datastats)./mean(groupls(4).datastats)*100; group(2).legend=['Incongruent'];
group(3).datastats = scaling*(groupos(1).datastats - groupos(2).datastats)./mean(groupos(2).datastats)*100; group(3).legend=['Congruent'];
group(4).datastats = scaling*(groupos(3).datastats - groupos(4).datastats)./mean(groupos(4).datastats)*100; group(4).legend=['Incongruent'];
i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group],opts_PSC);


%%
% Plot figures DELAY
opts_PSC.hmask = blkdiag(ones(2),ones(2),ones(2),ones(2));
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupld(cind)],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupld(iind)],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupod(cind)],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupod(iind)],[],opts_PM3Dcs,[]);
i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([groupld groupod],opts_PSC);


% Plot diffs DELAY
opts_PSC.hmask = blkdiag(ones(2),ones(2));
clear group; group(1:4) = Grp; group(1).freqband_stats = freqband_stats;
group(1).datastats = scaling*(groupld(1).datastats - groupld(2).datastats)./mean(groupld(2).datastats)*100; group(1).legend=['Congruent'];
group(2).datastats = scaling*(groupld(3).datastats - groupld(4).datastats)./mean(groupld(4).datastats)*100; group(2).legend=['Incongruent'];
group(3).datastats = scaling*(groupod(1).datastats - groupod(2).datastats)./mean(groupod(2).datastats)*100; group(3).legend=['Congruent'];
group(4).datastats = scaling*(groupod(3).datastats - groupod(4).datastats)./mean(groupod(4).datastats)*100; group(4).legend=['Incongruent'];
i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group],opts_PSC);




% Save figures
for i=7:12
%     figure(i);
%     set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})
end

%% Fig 3B2 - RT2 - All - PermuteSpect

clearvars -except wrkspc_buffer fv
vars_pull(fv);
currfigname = 'Fig3B2';

opts_PM3Dcs.do_sgolay = 1;

ffc_sfc_psd_mode = 3;
switch ffc_sfc_psd_mode
    case 1
        sfc_mode=22.4313101;    % RT
        perm_mode=22.4313111;
        script_name = 'plots_ffc_vs_Sch';
    case 2
        sfc_mode=2.4315100;
        perm_mode=2.4315110;
        script_name = 'plots_preferred';
    case 3
        sfc_mode=41.6313101;
        perm_mode=41.6313111;       % Incomplete
        script_name = 'plots_CFC';
end

[~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    

exclude_clipping = 0;
exclude_60 = 0;

perm2pls = 1;               % Instead of showing raw SFC values, show % cells successfully passing permutation test
        perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
        perm2pls_dophi = 0;

plot_SFC_spectra=0;

% Choose animals & stage
excludeL = 0; 
excludeO = 1; 
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
groupls = out.group;
clear out

% Choose animals & stage
excludeL = 0; 
excludeO = 1; 
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
groupld = out.group;
clear out

% Choose animals & stage
excludeL = 1; 
excludeO = 0; 
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
groupos = out.group;
clear out

% Choose animals & stage
excludeL = 1; 
excludeO = 0; 
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
groupod = out.group;
clear out

groupls(1).legend = 'Sample Congruent RT';
groupls(2).legend = 'Sample Incongruent RT';
groupld(1).legend = 'Delay Congruent RT';
groupld(2).legend = 'Delay Incongruent RT';
groupos(1).legend = 'Sample Congruent RT';
groupos(2).legend = 'Sample Incongruent RT';
groupod(1).legend = 'Delay Congruent RT';
groupod(2).legend = 'Delay Incongruent RT';

i=0;
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupls groupld],[],opts_PM3Dcs,[]);
% i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupos groupod],[],opts_PM3Dcs,[]);
% i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[],[],opts_PM3Dcs,[]);


% Save figures
% for i=1:2; figure(i); set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i}); end


%% Fig 4A - Boundary - Spectra


clearvars -except wrkspc_buffer fv
vars_pull(fv); %close all
currfigname = 'Fig4A';

% plotmode = 9;
% remove_dependent = 1;
gr_submode = 1:4;

ffc_sfc_psd_mode = 1;
switch ffc_sfc_psd_mode
    case 1
        sfc_mode=22.4413101;    % RT
        perm_mode=22.4413111;
        script_name = 'plots_ffc_vs_Sch';
        N_criteria=3;
        %freqband_stats = fv.freqband_stats_alpha;
        freqband_stats = fv.freqband_stats_lbeta;
        
    case 2
        sfc_mode=2.2415101;
        %perm_mode=2.2415111;
        script_name = 'plots_preferred';
        N_criteria=5;   % For Roy
        freqband_stats = fv.freqband_stats_hbeta;
    case 3
        sfc_mode=41.6413101;
        perm_mode=41.6413111;
        script_name = 'plots_CFC';
        N_criteria=3;
        freqband_stats = fv.freqband_stats_alpha;
        %freqband_stats = fv.freqband_stats_lbeta;
end

[~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);

    
group_override = get_group_overrides(fname_suffix_sfc, N_criteria, gr_submode);
if ~isempty(group_override); group = group_override; end


%opts_PSC.hmask = blkdiag(ones(2),ones(2),ones(2),ones(2));

% % Monkey L
excludeL = 0; 
excludeO = 1; 

exclude_60 = 0;

curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func(script_name);
groupls=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out

curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func(script_name);
groupld=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out

% % Monkey O
excludeL = 1; 
excludeO = 0; 

curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func(script_name);
groupos=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out

curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
    
push_workspace
out = script2func(script_name);
groupod=out.group;
wrkspc_buffer = out.wrkspc_buffer; clear out

cind = 1:2;
iind = 3:4;

switch ffc_sfc_psd_mode;
    case 1; units= 'FFC'; scaling=1;
    case 2; units= 'SFC'; scaling=1;
    case 3; units= 'PSD'; scaling=1;
end
%%

sample_delay_switch = 1;    % 1-Sample, 2-Delay

switch sample_delay_switch
    case 1
        grouplc = groupls;
        groupoc = groupos;
    case 2
        grouplc = groupld;
        groupoc = groupod;
end

if 1
    % Plot figures SAMPLE
    opts_PSC.hmask = blkdiag(true(2),true(2),true(2),true(2));
    i=0;
    % 50% morphs
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[grouplc([5,7,6])],[],opts_PM3Dcs,[]);
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupoc([5,7,6])],[],opts_PM3Dcs,[]);
    i=i+1; fign2; [h1, h, p, out.PSC] = plot_stats_custstruct([grouplc([5,6,7,8]) groupoc([5,6,7,8])],opts_PSC);

    % 60% morphs
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[grouplc([1,3,2])],[],opts_PM3Dcs,[]);
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupoc([1,3,2])],[],opts_PM3Dcs,[]);
    i=i+1; fign2; [h1, h, p, out.PSC] = plot_stats_custstruct([grouplc([1,2,3,4]) groupoc([1,2,3,4])],opts_PSC);

end
%%
do_percent = 1;

if do_percent
    % % Plot Diffs 50% Morphs 
    opts_PSC.hmask = blkdiag(true(2),true(2));
    clear group; group(1:4) = Grp; group(1).freqband_stats = freqband_stats;
    group(1).datastats = scaling*(grouplc(5).datastats - grouplc(6).datastats)./mean(grouplc(6).datastats)*100; group(1).legend=['50% rel.'];
    group(2).datastats = scaling*(grouplc(7).datastats - grouplc(8).datastats)./mean(grouplc(8).datastats)*100; group(2).legend=['50% irr.'];
    group(3).datastats = scaling*(groupoc(5).datastats - groupoc(6).datastats)./mean(groupoc(6).datastats)*100; group(3).legend=['50% rel.'];
    group(4).datastats = scaling*(groupoc(7).datastats - groupoc(8).datastats)./mean(groupoc(8).datastats)*100; group(4).legend=['50% irr.'];
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group],opts_PSC);

    % % Plot Diffs 60% Morphs 
    opts_PSC.hmask = blkdiag(true(2),true(2));
    clear group; group(1:4) = Grp; group(1).freqband_stats = freqband_stats;
    group(1).datastats = scaling*(grouplc(1).datastats - grouplc(2).datastats)./mean(grouplc(2).datastats)*100; group(1).legend=['60% rel.'];
    group(2).datastats = scaling*(grouplc(3).datastats - grouplc(4).datastats)./mean(grouplc(4).datastats)*100; group(2).legend=['60% irr.'];
    group(3).datastats = scaling*(groupoc(1).datastats - groupoc(2).datastats)./mean(groupoc(2).datastats)*100; group(3).legend=['60% rel.'];
    group(4).datastats = scaling*(groupoc(3).datastats - groupoc(4).datastats)./mean(groupoc(4).datastats)*100; group(4).legend=['60% irr.'];
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group],opts_PSC);
else
    % % Plot Diffs 50% Morphs 
    opts_PSC.hmask = blkdiag(true(2),true(2));
    clear group; group(1:4) = Grp; group(1).freqband_stats = freqband_stats;
    group(1).datastats = scaling*(grouplc(5).datastats - grouplc(6).datastats); group(1).legend=['50% rel.'];
    group(2).datastats = scaling*(grouplc(7).datastats - grouplc(8).datastats); group(2).legend=['50% irr.'];
    group(3).datastats = scaling*(groupoc(5).datastats - groupoc(6).datastats); group(3).legend=['50% rel.'];
    group(4).datastats = scaling*(groupoc(7).datastats - groupoc(8).datastats); group(4).legend=['50% irr.'];
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group],opts_PSC);

    % % Plot Diffs 60% Morphs 
    opts_PSC.hmask = blkdiag(true(2),true(2));
    clear group; group(1:4) = Grp; group(1).freqband_stats = freqband_stats;
    group(1).datastats = scaling*(grouplc(1).datastats - grouplc(2).datastats); group(1).legend=['60% rel.'];
    group(2).datastats = scaling*(grouplc(3).datastats - grouplc(4).datastats); group(2).legend=['60% irr.'];
    group(3).datastats = scaling*(groupoc(1).datastats - groupoc(2).datastats); group(3).legend=['60% rel.'];
    group(4).datastats = scaling*(groupoc(3).datastats - groupoc(4).datastats); group(4).legend=['60% irr.'];
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group],opts_PSC);
end

% % Save figures
% for i=[3,6]
%     figure(i);
%     set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})
% end

%% Fig 4C - Boundary - All - PermuteSpect

clearvars -except wrkspc_buffer fv
vars_pull(fv);
currfigname = 'Fig4C';

opts_PM3Dcs.do_sgolay = 0;
opts_PM3Dcs.stats_mode = 2;
inds = 1:4;

ffc_sfc_psd_mode = 1;
switch ffc_sfc_psd_mode
    case 1
        sfc_mode=22.4413101;    % RT
        perm_mode=22.4414111;
        script_name = 'plots_ffc_vs_Sch';
    case 2
        sfc_mode=2.4415100;     % Incomplete
        perm_mode=2.4415110;
        script_name = 'plots_preferred';
    case 3
        sfc_mode=41.6413101;
        perm_mode=41.6413111;       
        script_name = 'plots_CFC';
end

[~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    

exclude_60 = 0;

perm2pls = 1;               % Instead of showing raw SFC values, show % cells successfully passing permutation test
        perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
        perm2pls_dophi = 0;

plot_SFC_spectra=0;

% Choose animals & stage
excludeL = 0; 
excludeO = 0; 
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
groupls = out.group;
clear out

% Choose animals & stage
excludeL = 0; 
excludeO = 0; 
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
groupld = out.group;
clear out

% Choose animals & stage
excludeL = 1; 
excludeO = 0; 
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
groupos = out.group;
clear out

% Choose animals & stage
excludeL = 1; 
excludeO = 0; 
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
push_workspace; out = script2func(script_name); wrkspc_buffer = out.wrkspc_buffer;
groupod = out.group;
clear out


i=0;
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupls(inds) ],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupld(inds) ],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupos(inds) ],[],opts_PM3Dcs,[]);
i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupod(inds) ],[],opts_PM3Dcs,[]);

opts_PSC.hmask = blkdiag(true(2),true(2),true(2),true(2));
i=i+1; figure; [h1, h, p, out.PSC] = plot_stats_custstruct([groupls(inds) groupos(inds)],opts_PSC); % Sample stage
i=i+1; figure; [h1, h, p, out.PSC] = plot_stats_custstruct([groupld(inds) groupod(inds)],opts_PSC); % Delay stage

% Save figures
% for i=1:4; figure(i); set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i}); end


%% Huge regression (still gives shitty results)

% Globals (preserved after clear)
clearvars -except wrkspc_buffer fv

curr_stage=2;
    curr_stage_alt=3;

regress_pls_vs_sens = 0;

excludeL = 1;
excludeO = 0;

% Frequency ALPHA
do_alpha = 1;       % If 0 dose beta analysis; else alpha

% Sample, Sample
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
clear sensregstr; i=0;
i=i+1; sensregstr{i} = 'sqrt(sl1.*sr1)';
i=i+1; sensregstr{i} = 'sqrt(sl2.*sr2)';
i=i+1; sensregstr{i} = 'sqrt(sl9.*sr9)';
i=i+1; sensregstr{i} = 'adj11';         % Adjacency statistic (300-500Hz sfc)
i=i+1; sensregstr{i} = 'p11';
sensregnames = estimate_sens_names(sensregstr);

if ~excludeL && excludeO
    freqband_stats = [19.5 19.5];   % Sample Monkey L beta
elseif ~excludeO && excludeL
    freqband_stats = [17.5 17.5];        % Sample Monkey O beta
end
if do_alpha
    freqband_stats = [10.5 10.5];
end

push_workspace;
out = script2func('plots_preferred_unitpairs');
wrkspc_buffer = out.wrkspc_buffer;
out.sensregnames{4} = [out.sensregnames{4} ' sA'];
out.sensregnames{5} = ['FFC sA'];
sens_reg_sA = out.sens_reg;
sensregnames_sA = out.sensregnames; 
bad_regress_sA = out.bad_regress;
clear out


% Delay, Delay
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
clear sensregstr; i=0;
i=i+1; sensregstr{i} = 'sqrt(sla1.*sra1)';
i=i+1; sensregstr{i} = 'sqrt(sla2.*sra2)';
i=i+1; sensregstr{i} = 'sqrt(sla9.*sra9)';
i=i+1; sensregstr{i} = 'adj11';         % Adjacency statistic (300-500Hz sfc)
i=i+1; sensregstr{i} = 'p11';
sensregnames = estimate_sens_names(sensregstr);

freqband_stats = [16.6 16.6];
if do_alpha
    freqband_stats = [10.5 10.5];
end

push_workspace;
out = script2func('plots_preferred_unitpairs');
wrkspc_buffer = out.wrkspc_buffer;
out.sensregnames{4} = [out.sensregnames{4} ' dA'];
out.sensregnames{5} = ['FFC dA'];
sens_reg_dA = out.sens_reg;
sensregnames_dA = out.sensregnames; 
bad_regress_dA = out.bad_regress;
clear out


% Frequency BETA
do_alpha = 0;       % If 0 dose beta analysis; else alpha

% Sample, Sample
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
clear sensregstr; i=0;
i=i+1; sensregstr{i} = 'sqrt(sl1.*sr1)';
i=i+1; sensregstr{i} = 'sqrt(sl2.*sr2)';
i=i+1; sensregstr{i} = 'sqrt(sl9.*sr9)';
i=i+1; sensregstr{i} = 'adj11';         % Adjacency statistic (300-500Hz sfc)
i=i+1; sensregstr{i} = 'p11';
sensregnames = estimate_sens_names(sensregstr);

if ~excludeL && excludeO
    freqband_stats = [19.5 19.5];   % Sample Monkey L beta
elseif ~excludeO && excludeL
    freqband_stats = [17.5 17.5];        % Sample Monkey O beta
end
if do_alpha
    freqband_stats = [10.5 10.5];
end

push_workspace;
out = script2func('plots_preferred_unitpairs');
wrkspc_buffer = out.wrkspc_buffer;
out.sensregnames{4} = [out.sensregnames{4} ' sB'];
out.sensregnames{5} = ['FFC sB'];
sens_reg_sB = out.sens_reg;
sensregnames_sB = out.sensregnames; 
bad_regress_sB = out.bad_regress;
clear out


% Delay, Delay
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
clear sensregstr; i=0;
i=i+1; sensregstr{i} = 'sqrt(sla1.*sra1)';
i=i+1; sensregstr{i} = 'sqrt(sla2.*sra2)';
i=i+1; sensregstr{i} = 'sqrt(sla9.*sra9)';
i=i+1; sensregstr{i} = 'adj11';         % Adjacency statistic (300-500Hz sfc)
i=i+1; sensregstr{i} = 'p11';
sensregnames = estimate_sens_names(sensregstr);

freqband_stats = [16.6 16.6];
if do_alpha
    freqband_stats = [10.5 10.5];
end

push_workspace;
out = script2func('plots_preferred_unitpairs');
wrkspc_buffer = out.wrkspc_buffer;
out.sensregnames{4} = [out.sensregnames{4} ' dB'];
out.sensregnames{5} = ['FFC dB'];
sens_reg_dB = out.sens_reg;
sensregnames_dB = out.sensregnames; 
bad_regress_dB = out.bad_regress;
clear out


sens_all = [sens_reg_sA sens_reg_dA sens_reg_sB sens_reg_dB];
sens_names_all = [sensregnames_sA(:)' sensregnames_dA(:)' sensregnames_sB(:)' sensregnames_dB(:)'];
bad_all = [bad_regress_sA(:) bad_regress_dA(:) bad_regress_sB(:) bad_regress_dB(:) ];
bads = any(bad_all,2);

if 0
    % Plot all adj11 traces (justify removing the duplicates
    ind = find(cellfun(@(s) ~isempty(s),strfind(sens_names,'adj11')))
    figure; plot(sens_all(:,[4,9,14,19])); legend(sens_names{ind})
end

keep_ind = [1:3,6:8,4,9,5,10,15,20];    % Only keep unique data;
sens = sens_all(~bads,keep_ind);
snames = sens_names_all(keep_ind);


sens = zscore(sens);

[coeff, score, latent] = princomp(sens(:,9:10));
% figure; bar(latent);
% figure; bar(coeff(:,[1:2]))
% figure; plot(score(:,1),score(:,2),'.')
FFCA=score(:,1);

[coeff, score, latent] = princomp(sens(:,11:12));
% figure; bar(latent);
% figure; bar(coeff(:,[1:2]))
% figure; plot(score(:,1),score(:,2),'.')
FFCB=score(:,1);

[coeff, score, latent] = princomp(sens(:,7:8));
% figure; bar(latent);
% figure; bar(coeff(:,[1:2]))
% figure; plot(score(:,1),score(:,2),'.')
Adj=score(:,1);

[coeff, score, latent] = princomp(sens(:,[1:2,4:5]));
% figure; bar(latent);
% figure; bar(coeff(:,[1:2]))
% figure; plot(score(:,1),score(:,2),'.')
X=score(:,1);

[coeff, score, latent] = princomp(sens(:,[1:2]));
% figure; bar(latent);
% figure; bar(coeff(:,[1:2]))
% figure; plot(score(:,1),score(:,2),'.')
XA=score(:,1);

[coeff, score, latent] = princomp(sens(:,[4:5]));
% figure; bar(latent);
% figure; bar(coeff(:,[1:2]))
% figure; plot(score(:,1),score(:,2),'.')
XB=score(:,1);

[coeff, score, latent] = princomp(sens(:,[3,6]));
% figure; bar(latent);
% figure; bar(coeff(:,[1:2]))
% figure; plot(score(:,1),score(:,2),'.')
FR=score(:,1);

plots_regress([X,FR,Adj,FFCA],{'X','FR','Adj','FFCA'},[],'.');
plots_regress([X,FR,Adj,FFCB],{'X','FR','Adj','FFCA'},[],'.');

plots_regress([XA,XB,FR,Adj,FFCA],{'XA','XB','FR','Adj','FFCA'},[],'.');
plots_regress([XA,XB,FR,Adj,FFCB],{'XA','XB','FR','Adj','FFCA'},[],'.');

plots_regress([FFCA,Adj,FFCB],{'FFCA','adj','FFCB'},[],'.');



%% Figure 5A SFC Spectra vs cell groups

clearvars -except wrkspc_buffer fv
vars_pull(fv); 
currfigname = 'Fig5A';

sfc_mode=2.2015100;
perm_mode = 2.201511;
boot_mode = 2.201512;
    [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
    
    
animals_mode = 2;
if animals_mode == 1
    excludeL = 0; excludeO = 1; % Animals L
elseif animals_mode == 2
    excludeL = 1; excludeO = 0; % Animals O
end

plot_SFC_spectra = 0;
plot_SFC_spectra_stats = 0;

% % Group Setup
symmetry_mode=3;
[gr,mask] = group_setup_FR2D;       % Generate grouping structure
[gr.all_cells, gr.all_sens, gr.all_sens4, gr.single, gr.single_sens] = group_setup_other;
gr.rmr = score2grouping(gr.all,[gr.all.diff]);

group = gr.rmr;
group(1).legend='Multitasking';
group(2).legend='Specialist';
group(3).legend='Not sensitive';

% % % Uncomment this for deciders vs non-deciders
% clear group
% group(1) = group_merge(gr.rmr(1),gr.rmr(2)); group(1).legend = 'Any preference';
% group(2) = gr.rmr(3); group(2).legend = 'No preference';



% % Stages
profile on
curr_stage=2;
curr_stage_alt=3;
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
    curr_stage_stats = curr_stage_sfc;
push_workspace; out = script2func('plots_preferred'); wrkspc_buffer = out.wrkspc_buffer; 
groupss=out.group; clear out
profile off

% % Stages
curr_stage=2;
curr_stage_alt=3;
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
    curr_stage_stats = curr_stage_sfc;
push_workspace; out = script2func('plots_preferred'); wrkspc_buffer = out.wrkspc_buffer; 
groupsd=out.group; clear out

% % Stages
curr_stage=3;
curr_stage_alt=3;
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
    curr_stage_stats = curr_stage_sfc;
push_workspace; out = script2func('plots_preferred'); wrkspc_buffer = out.wrkspc_buffer; 
groupdd=out.group; clear out


append_ss = @(s) ['SS ',s];
append_sd = @(s) ['SD ',s];
append_dd = @(s) ['DD ',s];
groupss=myfunc(groupss,'legend',append_ss);
groupsd=myfunc(groupsd,'legend',append_sd);
groupdd=myfunc(groupdd,'legend',append_dd);

i=0;
if animals_mode == 1
%     i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupsd],[],opts_PM3Dcs,['bgrbgr']);
%     i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupss groupdd],[],opts_PM3Dcs,['cmkbgr']);
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupss ],[],opts_PM3Dcs,['bgrbgr']);
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupsd],[],opts_PM3Dcs,['bgrbgr']);
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupdd],[],opts_PM3Dcs,['bgrbgr']);
elseif animals_mode == 2
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupss ],[],opts_PM3Dcs,['bgrbgr']);
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupsd],[],opts_PM3Dcs,['bgrbgr']);
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupdd],[],opts_PM3Dcs,['bgrbgr']);
    
end

% Plot bargraph
% opts_PSC.hmask = blkdiag([1 0 1; 0 1 1; 1 1 1]);
% i=i+1; fign; [h1,hstats,pstats,out] = plot_stats_custstruct([ groupsd  ],opts_PSC);
% i=i+1; fign; [h1,hstats,pstats,out] = plot_stats_custstruct([ groupsd  ]);

% Plot bargraph 2
opts_PSC.hmask = blkdiag(true(3),true(3),true(3));
i=i+1; figure; [h1,hstats,pstats,out] = plot_stats_custstruct([groupss, groupsd, groupdd],opts_PSC);

%% Figure 5B SFC Spectra Regression

clearvars -except wrkspc_buffer fv
vars_pull(fv); 
currfigname = 'Fig5B';

sfc_mode=2.2015100;
perm_mode = 2.201511;
boot_mode = 2.201512;
    [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);

animals_mode = 1;
if animals_mode == 1
    excludeL = 0; excludeO = 1; % Animals L
elseif animals_mode == 2
    excludeL = 1; excludeO = 0; % Animals O
end

    

regress_pls_vs_sens = 1;

merge_schAB = 0;
i=0;
switch merge_schAB
    case 0
        i=i+1;sensregnames{i} = 'SchA';
        i=i+1;sensregnames{i} = 'SchB';
        ind = [2,3,4];
    
    case {1,2}
        i=i+1;sensregnames{i} = 'cPEV';
        ind = [2];
end
i=i+1;sensregnames{i} = 'FR';
i=i+1;sensregnames{i} = 'SFC';


% % Group Setup
symmetry_mode=3;
[gr,mask] = group_setup_FR2D;       % Generate grouping structure
[gr.all_cells, gr.all_sens, gr.all_sens4, gr.single, gr.single_sens] = group_setup_other;
gr.rmr = score2grouping(gr.all,[gr.all.diff]);
group = gr.rmr;


% % Stage 2/2
curr_stage=2;
curr_stage_alt=3;
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
    curr_stage_stats = curr_stage_sfc;
clear sensregstr; i=0;
switch merge_schAB
    case 0
    i=i+1; sensregstr{i} = 's1';
    i=i+1; sensregstr{i} = 's2';
%         i=i+1; sensregstr{i} = 'sqrt(s1.*s2)';
%         i=i+1; sensregstr{i} = 'sqrt(s1.^2+s2.^2)';
    case 1
    i=i+1; sensregstr{i} = 'sqrt(s1.*s2)';
    case 2
    i=i+1; sensregstr{i} = 's1+s2';
%     i=i+1; sensregstr{i} = 'sqrt(s1+s2)';
end
i=i+1; sensregstr{i} = 's9';    % Alternate firing rate (for delay stage, most likely - we should use this)
i=i+1; sensregstr{i} = 'p11';
push_workspace; out = script2func('plots_preferred'); wrkspc_buffer = out.wrkspc_buffer; 
% fign; plot_reg_custstruct([out.reg_struct])
regss = out.reg_struct; clear out


% % Stage 2/3
curr_stage=2;
curr_stage_alt=3;
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
    curr_stage_stats = curr_stage_sfc;
clear sensregstr; i=0;
switch merge_schAB
    case 0
    i=i+1; sensregstr{i} = 's1';
    i=i+1; sensregstr{i} = 's2';
%         i=i+1; sensregstr{i} = 'sqrt(s1.*s2)';
%         i=i+1; sensregstr{i} = 'sqrt(s1.^2+s2.^2)';
    case 1
    i=i+1; sensregstr{i} = 'sqrt(s1.*s2)';
    case 2
    i=i+1; sensregstr{i} = 's1+s2';
%     i=i+1; sensregstr{i} = 'sqrt(s1+s2)';
end
i=i+1; sensregstr{i} = 'sa9';    % Alternate firing rate (for delay stage, most likely - we should use this)
i=i+1; sensregstr{i} = 'p11';
push_workspace; out = script2func('plots_preferred'); wrkspc_buffer = out.wrkspc_buffer; 
% fign; plot_reg_custstruct([out.reg_struct])
regsd = out.reg_struct; clear out


% Stage 3/3
clear sensregstr; i=0;
switch merge_schAB
    case 0
    i=i+1; sensregstr{i} = 'sa1';
    i=i+1; sensregstr{i} = 'sa2';
%         i=i+1; sensregstr{i} = 'sqrt(sa1.*sa2)';
%         i=i+1; sensregstr{i} = 'sqrt(sa1.^2+sa2.^2)';
    case 1
    i=i+1; sensregstr{i} = 'sqrt(sa1.*sa2)';
    case 2
    i=i+1; sensregstr{i} = 'sa1+sa2';
%     i=i+1; sensregstr{i} = 'sqrt(sa1+sa2)';
end
i=i+1; sensregstr{i} = 'sa9';    % Alternate firing rate (for delay stage, most likely - we should use this)
i=i+1; sensregstr{i} = 'p11';
push_workspace; out = script2func('plots_preferred'); wrkspc_buffer = out.wrkspc_buffer; 
% fign; plot_reg_custstruct([out.reg_struct])
regdd = out.reg_struct; clear out

%figure; plot_reg_custstruct([regss regsd regdd],opts_PSC)


temp=strcat({regss.legend_arr},' SS'); [regss.legend_arr] = temp{:};
temp=strcat({regsd.legend_arr},' SD'); [regsd.legend_arr] = temp{:};
temp=strcat({regdd.legend_arr},' DD'); [regdd.legend_arr] = temp{:};
fign; plot_reg_custstruct([regss(ind) regsd(ind) regdd(ind)]);


% 
% i=0;
% i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupls groupld],[],opts_PM3Dcs,['bgbg']); set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})
% i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupos groupod],[],opts_PM3Dcs,['bgbg']); set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})



%% Figure 5C FFC Spectra Regression

clearvars -except wrkspc_buffer fv
vars_pull(fv); 
currfigname = 'Fig5C';

sfc_mode=22.4013100;
perm_mode=22.401311;
boot_mode = 22.401312;
sfc_mode_forreal = 2.2015100;  % We just need this to get spikerate!
    [~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);

animals_mode = 2;
if animals_mode == 1
    excludeL = 0; excludeO = 1; % Animals L
elseif animals_mode == 2
    excludeL = 1; excludeO = 0; % Animals O
end

do_alpha = 1;
%exclude_60 = 0;

% Sample Sample
curr_stage_sfc=2;
    curr_stage_badclipping = curr_stage_sfc;
clear sensregstr; i=0;
i=i+1; sensregstr{i} = 'sqrt(sl1.*sr1)';
i=i+1; sensregstr{i} = 'sqrt(sl2.*sr2)';
i=i+1; sensregstr{i} = 'sqrt(sl9.*sr9)';
i=i+1; sensregstr{i} = 'adj11';         % Adjacency statistic (300-500Hz sfc)
i=i+1; sensregstr{i} = 'p11';
sensregnames = estimate_sens_names(sensregstr);

if ~excludeL && excludeO
    freqband_stats = [19.5 19.5];   % Sample Monkey L beta
elseif ~excludeO && excludeL
    freqband_stats = [17.5 17.5];        % Sample Monkey O beta
end
if do_alpha
    freqband_stats = [10.5 10.5];
end

regress_pls_vs_sens = 1;

push_workspace;
out = script2func('plots_preferred_unitpairs');
wrkspc_buffer = out.wrkspc_buffer;
regss = out.reg_struct; clear out

% Sample Delay
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
clear sensregstr; i=0;
i=i+1; sensregstr{i} = 'sqrt(sl1.*sr1)';
i=i+1; sensregstr{i} = 'sqrt(sl2.*sr2)';
i=i+1; sensregstr{i} = 'sqrt(sla9.*sra9)';
i=i+1; sensregstr{i} = 'adj11';         % Adjacency statistic (300-500Hz sfc)
i=i+1; sensregstr{i} = 'p11';
sensregnames = estimate_sens_names(sensregstr);

freqband_stats = [16.6 16.6];
if do_alpha
    freqband_stats = [10.5 10.5];
end

push_workspace;
out = script2func('plots_preferred_unitpairs');
wrkspc_buffer = out.wrkspc_buffer;
regsd = out.reg_struct; clear out

% Delay Delay
curr_stage_sfc=3;
    curr_stage_badclipping = curr_stage_sfc;
clear sensregstr; i=0;
i=i+1; sensregstr{i} = 'sqrt(sla1.*sra1)';
i=i+1; sensregstr{i} = 'sqrt(sla2.*sra2)';
i=i+1; sensregstr{i} = 'sqrt(sla9.*sra9)';
i=i+1; sensregstr{i} = 'adj11';         % Adjacency statistic (300-500Hz sfc)
i=i+1; sensregstr{i} = 'p11';
sensregnames = estimate_sens_names(sensregstr);

freqband_stats = [16.6 16.6];
if do_alpha
    freqband_stats = [10.5 10.5];
end

push_workspace;
out = script2func('plots_preferred_unitpairs');
wrkspc_buffer = out.wrkspc_buffer;
regdd = out.reg_struct; clear out




ind=2:3

temp=strcat({regss.legend_arr},' SS'); [regss.legend_arr] = temp{:};
temp=strcat({regsd.legend_arr},' SD'); [regsd.legend_arr] = temp{:};
temp=strcat({regdd.legend_arr},' DD'); [regdd.legend_arr] = temp{:};
fign; plot_reg_custstruct([regss(ind) regsd(ind) regdd(ind)]);


% 
% i=0;
% i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupls groupld],[],opts_PM3Dcs,['bgbg']); set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})
% i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[groupos groupod],[],opts_PM3Dcs,['bgbg']); set(gcf,'PaperPositionMode','auto'); print(gcf,'-dpng','-r200',savenames{i})



