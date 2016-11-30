
% New plotting function that uses generalized code
%   See plots_generalized for template or funcs_generalized for useful
%   functions

%%
error('This is meant to run in cell mode. Dont F5!');


%%
% % For loop for saving figs
if ~exist('currfname'); currfname = 'figs_batch_2_0_ctgs'; end
if ~exist('currfigname'); currfigname = 'Fig6'; end
savenames={'fig1','fig2','fig3','fig4','fig5','fig6','fig7','fig8','fig9','fig10','fig11','fig12','fig13','fig14','fig15','fig16','fig17','fig18','fig19','fig20','fig21','fig22','fig23','fig24'};
mydate = datestr(datenum(date),'yy/mm/dd'); mydate = strrep(mydate,'/','');
c=clock;
sp = ['d' mydate '_t' num2str(c(4),'%10.2d') '' num2str(c(5),'%10.2d') '' num2str(round(c(6)),'%10.2d')];
sp = [sp '__' currfname '_' currfigname];
basepath = '.';
% basepath = '~/figs_tosave';
mkdir(fullfile(basepath,sp));
do_pdf = 0;
for i=[1:4]
    figure(i); %ylim([0 0.175])
    title('');
    ylabel('');
    %set(gcf,'Position',[ 421    60   582   746]);
    %xlim([-1.5 2.2]);
    %ylabel('Avg z-score |\Delta FFC|')
    if ~do_pdf
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-dpng','-r200',fullfile(basepath,sp,savenames{i}))
        %print(gcf,'-dpdf',fullfile(basepath,sp,savenames{i}))
        %print(gcf,'-dpng',fullfile('.',sp,savenames{i}))
    else
%         set(gcf,'PaperPositionMode','auto');
        saveas(gcf,fullfile(basepath,sp,savenames{i}),'pdf');
        %print(gcf,'-dpdf',fullfile(basepath,sp,savenames{i}))
    end
    %close
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

% Clear
clearvars -except wrkspc_buffer fv fv3

addpath('figs_batch_3_0_supporting');
addpath(genpath('funcs_plots_preferred'));
addpath(genpath('funcs_generalized'));

% Load defaults
fv3 = figs_batch_defaults_3; 

% Create wrkspc_buffer structure if not already present
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end


%% Preload all SFC data

if 1
    warning('This will load all data might fill up memory. Press any key to continue.');
    
    
    % Setup variables
    %sfc_mode = [22.4413101, 22.4513101, 22.4914103];
    sfc_mode = [22.4413101, 22.4513101, 22.4914103];
    perm_mode = [22.4014111,22.4015111,22.4017111,...
        22.40141111,22.40151111,22.40171111,...
        22.4514111,22.4515111,...
        22.45141111,22.45151111,22.45171111,...
        ];
    path_matfiles_store = getpath('path_matfiles_store');


% 
%     % Load the pairs data (coherence)
%     for i = 1:length(sfc_mode)
%         % SFC
%         buff_fieldname = ['sfc_' mode2modename(sfc_mode(i))];
%         buff_func = @() func_recalc_pairs(stage_range,file_range,sfc_mode(i));
%         wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,0,0);
%     end

    for i = 1:length(perm_mode)
        % Load perm_mode tests
        buff_mode = perm_mode(i);
        buff_fieldname = ['per_' mode2modename(buff_mode)];
        buff_func = @() func_recalc_perms(buff_mode,[2,3],1:79);
        wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,0,0);

    end


    clear sfc_mode perm_mode

end

%% Fig 6c - Ctgs1v2 - All - PermuteSpect - Scatterplot

clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig6c';

do_both_animals=1;
script_name = 'plots_generalized';


% SFC mode
sfc_mode = 22.4014111;
perm_mode = sfc_mode;
% sfc_mode = 41.6013111;
% perm_mode = 41.6013111;


    % Pls switches
    freqband_stats = [16 20];

    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    % More pls switches
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.do_diff = 1;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
    % Permutation test options
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;


    % Create group template
    mycrit = [2*ones(1,5)];
    grt = Grp;
    grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
    clear grt2
    i=0;
    i=i+1; grt2(i)=grt; grt2(i).criteria(1:2)=[1 0]; grt2(i).ctgs=1; grt2(i).legend = 'Pref SchA';
    i=i+1; grt2(i)=grt; grt2(i).criteria(1:2)=[0 1]; grt2(i).ctgs=1; grt2(i).legend = 'Pref SchB';
    i=i+1; grt2(i)=grt; grt2(i).criteria(5)=[1]; grt2(i).ctgs=1; grt2(i).legend = 'Rule';
    i=i+1; grt2(i)=grt; grt2(i).criteria(1:2)=[1 1]; grt2(i).ctgs=1; grt2(i).legend = 'Multitasker';
    %i=i+1; gr(i)=grt; gr(i).criteria(1:2)=[0 0]; gr(i).ctgs=1;   % Ctg1-2 non-deciders
    gr = grt2;

    % Create group template for plus-minus counting
    mycrit = [2*ones(1,10)];
    grt = Grp;
    grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
    clear grt2
    i=0;
    i=i+1; grt2(i)=grt; grt2(i).criteria(1:4)=[1 0 0 0]; grt2(i).ctgs=1; grt2(i).legend = 'Pref SchA+';
    i=i+1; grt2(i)=grt; grt2(i).criteria(1:4)=[0 1 0 0]; grt2(i).ctgs=1; grt2(i).legend = 'Pref SchA-';
    i=i+1; grt2(i)=grt; grt2(i).criteria(1:4)=[0 0 1 0]; grt2(i).ctgs=1; grt2(i).legend = 'Pref SchB+';
    i=i+1; grt2(i)=grt; grt2(i).criteria(1:4)=[0 0 0 1]; grt2(i).ctgs=1; grt2(i).legend = 'Pref SchB-';
    i=i+1; grt2(i)=grt; grt2(i).criteria(9:10)=[1 0]; grt2(i).ctgs=1; grt2(i).legend = 'Rule+';
    i=i+1; grt2(i)=grt; grt2(i).criteria(9:10)=[0 1]; grt2(i).ctgs=1; grt2(i).legend = 'Rule-';
    gr_plusminus = grt2;

    
    % % % Animal L
    % Stage 2
    stage_sfc=2; stage_perm=2;
    opts_exclude.excludeL = 0;
    opts_exclude.excludeO = 1; 
    [wrkspc_buffer, group, pls, pls_stats, bad_any] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,gr,opts_exclude,opts_pls,opts_perm);
    opts_perm.split_plusminus=1; [~, group_pm, ~, ~, ~] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,gr_plusminus,opts_exclude,opts_pls,opts_perm); opts_perm.split_plusminus=0; fprintf(['Sample_L - N SchA +ve:' num2str(group_pm(1).numcells) ' N SchA -ve:' num2str(group_pm(2).numcells) ' N SchB +ve:' num2str(group_pm(3).numcells) ' N SchB -ve:' num2str(group_pm(4).numcells) ' N CueA:' num2str(group_pm(5).numcells) ' N CueB:' num2str(group_pm(6).numcells) '\n']); 
    figure; x = pls_stats(~bad_any,1); y = pls_stats(~bad_any,2); plot(x,y,'k.');
    hold on; plot_scattergroups(x,y,group,~bad_any);
    xlabel('FFC_{Cat}-FFC_{Dog}'); ylabel('FFC_{Tad}-FFC_{Goc}'); myx = xlim; myy = ylim; hold on; plot([-1 1],[0,0],'k','LineWidth',2); plot([0,0],[-1 1],'k','LineWidth',2); xlim(myx); ylim(myy);
    
    %
    % Stage 3
    stage_sfc=3; stage_perm=3;
    opts_exclude.excludeL = 0;
    opts_exclude.excludeO = 1; 
    [wrkspc_buffer, group, pls, pls_stats, bad_any] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,gr,opts_exclude,opts_pls,opts_perm);
    opts_perm.split_plusminus=1; [~, group_pm, ~, ~, ~] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,gr_plusminus,opts_exclude,opts_pls,opts_perm); opts_perm.split_plusminus=0; fprintf(['Delay_L - N SchA +ve:' num2str(group_pm(1).numcells) ' N SchA -ve:' num2str(group_pm(2).numcells) ' N SchB +ve:' num2str(group_pm(3).numcells) ' N SchB -ve:' num2str(group_pm(4).numcells) ' N CueA:' num2str(group_pm(5).numcells) ' N CueB:' num2str(group_pm(6).numcells) '\n']); 
    figure; x = pls_stats(~bad_any,1); y = pls_stats(~bad_any,2); plot(x,y,'k.');
    hold on; plot_scattergroups(x,y,group,~bad_any);
    xlabel('FFC_{Cat}-FFC_{Dog}'); ylabel('FFC_{Tad}-FFC_{Goc}'); myx = xlim; myy = ylim; hold on; plot([-1 1],[0,0],'k','LineWidth',2); plot([0,0],[-1 1],'k','LineWidth',2); xlim(myx); ylim(myy);
    
    %
    % % % Animal O
    % Stage 2
    stage_sfc=2; stage_perm=2;
    opts_exclude.excludeL = 1;
    opts_exclude.excludeO = 0; 
    [wrkspc_buffer, group, pls, pls_stats, bad_any] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,gr,opts_exclude,opts_pls,opts_perm);
    opts_perm.split_plusminus=1; [~, group_pm, ~, ~, ~] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,gr_plusminus,opts_exclude,opts_pls,opts_perm); opts_perm.split_plusminus=0; fprintf(['Sample_O - N SchA +ve:' num2str(group_pm(1).numcells) ' N SchA -ve:' num2str(group_pm(2).numcells) ' N SchB +ve:' num2str(group_pm(3).numcells) ' N SchB -ve:' num2str(group_pm(4).numcells) ' N CueA:' num2str(group_pm(5).numcells) ' N CueB:' num2str(group_pm(6).numcells) '\n']); 
    figure; x = pls_stats(~bad_any,1); y = pls_stats(~bad_any,2); plot(x,y,'k.');
    hold on; plot_scattergroups(x,y,group,~bad_any);
    xlabel('FFC_{Cat}-FFC_{Dog}'); ylabel('FFC_{Tad}-FFC_{Goc}'); myx = xlim; myy = ylim; hold on; plot([-1 1],[0,0],'k','LineWidth',2); plot([0,0],[-1 1],'k','LineWidth',2); xlim(myx); ylim(myy);
    
    %
    % Stage 3
    stage_sfc=3; stage_perm=3;
    opts_exclude.excludeL = 1;
    opts_exclude.excludeO = 0; 
    [wrkspc_buffer, group, pls, pls_stats, bad_any] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,gr,opts_exclude,opts_pls,opts_perm);
    opts_perm.split_plusminus=1; [~, group_pm, ~, ~, ~] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,gr_plusminus,opts_exclude,opts_pls,opts_perm); opts_perm.split_plusminus=0; fprintf(['Delay_O - N SchA +ve:' num2str(group_pm(1).numcells) ' N SchA -ve:' num2str(group_pm(2).numcells) ' N SchB +ve:' num2str(group_pm(3).numcells) ' N SchB -ve:' num2str(group_pm(4).numcells) ' N CueA:' num2str(group_pm(5).numcells) ' N CueB:' num2str(group_pm(6).numcells) '\n']); 
    figure; x = pls_stats(~bad_any,1); y = pls_stats(~bad_any,2); plot(x,y,'k.');
    hold on; plot_scattergroups(x,y,group,~bad_any);
    xlabel('FFC_{Cat}-FFC_{Dog}'); ylabel('FFC_{Tad}-FFC_{Goc}'); myx = xlim; myy = ylim; hold on; plot([-1 1],[0,0],'k','LineWidth',2); plot([0,0],[-1 1],'k','LineWidth',2); xlim(myx); ylim(myy);
    
    %%
    % % % Both
    % Stage 2
    stage_sfc=2; stage_perm=2;
    opts_exclude.excludeL = 0;
    opts_exclude.excludeO = 0; 
    [wrkspc_buffer, group, pls, pls_stats, bad_any] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,gr,opts_exclude,opts_pls,opts_perm);
    opts_perm.split_plusminus=1; [~, group_pm, ~, ~, ~] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,gr_plusminus,opts_exclude,opts_pls,opts_perm); opts_perm.split_plusminus=0; fprintf(['Sample_Both - N SchA +ve:' num2str(group_pm(1).numcells) ' N SchA -ve:' num2str(group_pm(2).numcells) ' N SchB +ve:' num2str(group_pm(3).numcells) ' N SchB -ve:' num2str(group_pm(4).numcells) ' N CueA:' num2str(group_pm(5).numcells) ' N CueB:' num2str(group_pm(6).numcells) '\n']); 
    figure; x = pls_stats(~bad_any,1); y = pls_stats(~bad_any,2); plot(x,y,'k.');
    hold on; plot_scattergroups(x,y,group,~bad_any);
    xlabel('FFC_{Cat}-FFC_{Dog}'); ylabel('FFC_{Tad}-FFC_{Goc}'); myx = xlim; myy = ylim; hold on; plot([-1 1],[0,0],'k','LineWidth',2); plot([0,0],[-1 1],'k','LineWidth',2); xlim(myx); ylim(myy);
    %
    % Stage 3
    stage_sfc=3; stage_perm=3;
    opts_exclude.excludeL = 0;
    opts_exclude.excludeO = 0; 
    [wrkspc_buffer, group, pls, pls_stats, bad_any] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,gr,opts_exclude,opts_pls,opts_perm);
    opts_perm.split_plusminus=1; [~, group_pm, ~, ~, ~] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,gr_plusminus,opts_exclude,opts_pls,opts_perm); opts_perm.split_plusminus=0; fprintf(['Delay_Both - N SchA +ve:' num2str(group_pm(1).numcells) ' N SchA -ve:' num2str(group_pm(2).numcells) ' N SchB +ve:' num2str(group_pm(3).numcells) ' N SchB -ve:' num2str(group_pm(4).numcells) ' N CueA:' num2str(group_pm(5).numcells) ' N CueB:' num2str(group_pm(6).numcells) '\n']); 
    figure; x = pls_stats(~bad_any,1); y = pls_stats(~bad_any,2); plot(x,y,'k.');
    hold on; plot_scattergroups(x,y,group,~bad_any);
    xlabel('FFC_{Cat}-FFC_{Dog}'); ylabel('FFC_{Tad}-FFC_{Goc}'); myx = xlim; myy = ylim; hold on; plot([-1 1],[0,0],'k','LineWidth',2); plot([0,0],[-1 1],'k','LineWidth',2); xlim(myx); ylim(myy);
    
%% Fig 6e2 - Ctgs - FFC/Units plot rel vs irrel plot
% Doesn't combine sample & delay, but does work with units
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig6e2';

do_units = 1;

% Plot switches
plot_on_spect = 1;
plot_on_bargraph = 1;



% SFC mode
if do_units == 0
    sfc_mode = 22.4014111;
    perm_mode = 22.4014111;
    % sfc_mode = 22.401511111;  % Partial
    % sfc_mode = 2.3015111;
    % perm_mode = 2.3015111;
    % perm_mode = 5;
    % sfc_mode = 41.6514111;
    % perm_mode = 41.6514111;
elseif do_units == 1
    sfc_mode = 52.7002010;
    perm_mode = 52.7000010;
    plot_on_bargraph = 0;
elseif do_units == 2
    sfc_mode = 52.7000010;
    perm_mode = 52.7000010;
    plot_on_spect = 0;
end


% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 3;

if do_units==1
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
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
    opts_pls.collapse_pls_to_days = 0;
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dcs.stats_mode = 0;
    opts_PM3Dcs.remove_dependent = 0;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    
    warning('To do - compare sch sensitive electrodes ability to attenuate irrelevant ctg response to other electrodes.');


        
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

% Build group

% Group based on perm_mode
group = group0;
sp = ~bad_any(:);

i=0;
i=i+1; group(i).legend = 'Cat vs. Dog rel.'
i=i+1; group(i).legend = 'Goc vs. Tad rel.'
i=i+1; group(i).legend = 'Cat vs. Dog irr.'
i=i+1; group(i).legend = 'Goc vs. Tad irr.'
i=i+1; group(i).legend = 'SchA vs. SchB'

% i=0;
% i=i+1; group(i).legend = 'C/D Rel.'
% i=i+1; group(i).legend = 'F/T Rel.'
% i=i+1; group(i).legend = 'C/D Irrel.'
% i=i+1; group(i).legend = 'F/T Irrel.'
% i=i+1; group(i).legend = 'Sch. A/B'

% Load data into groups
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);

% Test Plot groups, at last
if plot_on_spect
    opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
    
    
    inds = 1:length(group);
    inds = [1,3];
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    
    inds = [2,4];
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    
%     figure;
%     inds = 1:4;
%     [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D);
%     figure;
%     inds = 5:8;
%     [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D);
end



if plot_on_bargraph
    fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group([1,3,2,4]),opts_PSC);
end



%% Fig 6f - Plot histograms
% Doesn't combine sample & delay, but does work with units
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig6f';

do_units = 0;

% Plot switches
plot_on_spect = 0;
plot_on_bargraph = 0;
plot_on_hist = 1;

% SFC mode
if do_units == 0
    sfc_mode = 22.401811101;
    perm_mode = 22.401811101;
    % sfc_mode = 2.3015111;
    % perm_mode = 2.3015111;
    % perm_mode = 5;
    % sfc_mode = 41.6514111;
    % perm_mode = 41.6514111;
elseif do_units == 1
    sfc_mode = 52.700201001;
    perm_mode = 52.700001001;
    plot_on_bargraph = 0;
elseif do_units == 2
    sfc_mode = 52.700001001;
    perm_mode = 52.700001001;
    plot_on_spect = 0;
end


% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 3;

if do_units==1
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
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 1;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dcs.stats_mode = 0;
    opts_PM3Dcs.remove_dependent = 0;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    
        
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

% Build group
% Load sp's
[wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm,opts_exclude);
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
i=i+1; group(i)=grt; group(i).criteria(1:2)=[1 2]; group(i).ctgs=1;   % Ctg1-2 deciders
i=i+1; group(i)=grt; group(i).criteria(1:2)=[0 2]; group(i).ctgs=1;   % Ctg1-2 non-deciders
i=i+1; group(i)=grt; group(i).criteria(1:2)=[2 1]; group(i).ctgs=2;   % Ctg1-2 deciders
i=i+1; group(i)=grt; group(i).criteria(1:2)=[2 0]; group(i).ctgs=2;   % Ctg1-2 non-deciders
% i=i+1; group(i)=grt; group(i).criteria(5)=[1]; group(i).ctgs=5;   % Ctg1-2 deciders
% i=i+1; group(i)=grt; group(i).criteria(5)=[0]; group(i).ctgs=5;   % Ctg1-2 non-deciders


% Calculate legend entries
group = group.query_legend(group0);

% Load data into groups
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);

% Test Plot groups, at last
if plot_on_spect
    opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
    
    
    inds = 1:length(group);
%     inds = [1,3];
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    
%     figure;
%     inds = 1:4;
%     [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D);
%     figure;
%     inds = 5:8;
%     [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D);
end



if plot_on_bargraph
    fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group([1,3,2,4]),opts_PSC);
end

if plot_on_hist
    %% Plot histograms
    fontsize = 40;
    i=0;
    i=i+1; group(i).legend = 'Sig';
    i=i+1; group(i).legend = 'Not-Sig';
    i=i+1; group(i).legend = 'Sig';
    i=i+1; group(i).legend = 'Not-Sig';
%     i=i+1; group(i).legend = 'Sig';
%     i=i+1; group(i).legend = 'Not-Sig';


    hist_settings.nbins_or_edges = 50;
    [h1] = plot_histplot(group,'plot_type','hist_paired','hist_settings',hist_settings);      % Using new plot histogram code.
    
    
end


%% Fig 8c - Do specialist cells exhibit enhanced boundary response for their preferred boundary?
% Redone from figs_batch_2

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig8c';
plot_on_spect = 1;
plot_on_bargraph = 1;
do_scatterplot = 0;

% Set up default structures
    % Pls switches
    freqband_stats = [16 20];
    
    % Set up pls options
    % Unit exclusion
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    % More pls switches
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_dophi = 0;
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.do_diff = 1;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))

    % Perm opts for getting sp's
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    

% Bar plotting struct
opts_PM3Dcs.remove_dependent=1;
opts_PSC.remove_dependent=1;
    
% Setup group
myN_criteria = 5;
grt = Grp;
grt.criteria = [2*ones(1,myN_criteria)]; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
grt.xlims_desired = [0 120];


% Criteria
mycrit_deciderA = [1 0 2 2 2;];   % Sch A and/or B deciders
mycrit_deciderB = [0 1 2 2 2;];   % Sch A and/or B deciders
mycrit_nondecider = [0 0 2 2 2];

% Ctgs
myctgs_bndryA = [1];      % Sch A Bndry Rel
myctgs_bndryB = [2];      % Sch B Bndry Rel 
myctgs_bndryB_irr = [3];  % Sch B Bndry IrrRel (i.e. Sch A cued)
myctgs_bndryA_irr = [4];  % Sch A Bndry IrrRel (i.e. Sch B cued)


clear mygr
i=0;

% Sch A Deciders; SchA Bndry Rel vs SchB Bndry Rel vs SchB Bndry Irrel
i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryA,size(mycrit,1),1)]; mygr(i).legend = 'dAbA';
i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryB,size(mycrit,1),1)]; mygr(i).legend = 'dAbB';
i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryB_irr,size(mycrit,1),1)]; mygr(i).legend = 'dAbBirr';

% Sch B Deciders; SchB Bndry Rel vs SchA Bndry Rel vs SchA Bndry Irrel
i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryB,size(mycrit,1),1)]; mygr(i).legend = 'dBbB';
i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryA,size(mycrit,1),1)]; mygr(i).legend = 'dBbA';
i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryA_irr,size(mycrit,1),1)]; mygr(i).legend = 'dBbAirr';

% 
% Set up sim specific parameters

% Set stages and excludes
stage_sfc = 3;
opts_exclude.excludeL = 0; 
opts_exclude.excludeO = 0; 

% Load data and build group for Fig 8c
[wrkspc_buffer, group] = Fig_8c_Bndry_vs_Ctgs2(wrkspc_buffer,stage_sfc,freqband_stats,mygr,opts_exclude,opts_pls,opts_perm);

if plot_on_spect
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[group(1:3)],[],opts_PM3Dcs,[]);
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[group(4:6)],[],opts_PM3Dcs,[]);
end

if plot_on_bargraph
    opts_PSC.hmask = blkdiag(true(2),true(2));
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group([1,2,4,5]),opts_PSC);
end




% Do scatterplot
if do_scatterplot
    % Set stages and excludes
    do_abs = 0;

    % Load data and build group for Fig 8c
    [wrkspc_buffer, group] = Fig_8c_Bndry_vs_Ctgs2_scatter(wrkspc_buffer,stage_sfc,freqband_stats,mygr,opts_exclude,opts_pls,opts_perm,do_abs);

end

%% Fig 10 - Do cells increase their FFC in response to overall task relevance of FFC modulation.
% Redone from figs_batch_2

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig10';
plot_on_spect = 1;
plot_on_bargraph = 0;

% Set up default structures
    % Pls switches
    freqband_stats = [16 20];
    
    % Set up pls options
    % Unit exclusion
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    % More pls switches
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_dophi = 0;
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))

    % Perm opts for getting sp's
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    

% Bar plotting struct
opts_PM3Dcs.remove_dependent=1;
opts_PSC.remove_dependent=0;
    
% Setup group
myN_criteria = 5;
grt = Grp;
grt.criteria = [2*ones(1,myN_criteria)]; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 9;
grt.xlims_desired = [0 120];


% Criteria
mycrit_multi = [1 1 2 2 2;];   % Sch A and/or B deciders
mycrit_deciderA = [1 0 2 2 2;];   % Sch A and/or B deciders
mycrit_deciderB = [0 1 2 2 2;];   % Sch A and/or B deciders
mycrit_nondecider = [0 0 2 2 2];

clear mygr
i=0;

% Sch A Deciders; SchA Bndry Rel vs SchB Bndry Rel vs SchB Bndry Irrel
% i=i+1; mycrit = mycrit_multi; mygr(i) = grt; mygr(i).criteria=[mycrit]; 
i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; 
i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; 
i=i+1; mycrit = mycrit_nondecider; mygr(i) = grt; mygr(i).criteria=[mycrit]; 


% 
% Set up sim specific parameters
sfc_mode =  22.40131110;
perm_mode = 22.40141110;
% sfc_mode =  41.60131110;
% perm_mode = 41.60141110;


% Set stages and excludes
stage_perm = 3;
stage_sfc = 3;
opts_exclude.excludeL = 1; 
opts_exclude.excludeO = 0; 

% Load data and build group for Fig 8c
[wrkspc_buffer, group] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,mygr,opts_exclude,opts_pls,opts_perm);


if plot_on_spect
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[group(:)],[],opts_PM3Dcs,[]);
end

if plot_on_bargraph
    opts_PSC.hmask = blkdiag(true(2),true(2));
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group(:),opts_PSC);
end




%% Fig 11 - Ctg deciders Mag(C) vs Phi(C)
% Redone from figs_batch_2

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig11';
plot_on_spect = 1;
plot_on_bargraph = 1;
plot_on_scatter = 0;

% Set up default structures
    % Pls switches
    freqband_stats = [16 20];
    
    % Set up pls options
    % Unit exclusion
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    % More pls switches
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        %opts_pls.perm2pls_dophi = 1;
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))

    % Perm opts for getting sp's
    opts_perm.do_bh0 = 1;
    %opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    

% Bar plotting struct
opts_PM3Dcs.remove_dependent=1;
opts_PSC.remove_dependent=0;

opts_PM3Dcs.stats_mode = 2;
    
% Setup group
myN_criteria = 5;
grt = Grp;
grt.criteria = [2*ones(1,myN_criteria)]; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 10;
grt.xlims_desired = [0 200]; grt.interpreter = 'latex';


% Criteria
mycrit_multi = [1 1 2*ones(1,myN_criteria - 2);];   % Sch A and/or B deciders

% Criteria - exclusive
mycrit_deciderA = [1 0 2*ones(1,myN_criteria - 2);];   % Sch A and/or B deciders
mycrit_deciderB = [0 1 2*ones(1,myN_criteria - 2);];   % Sch A and/or B deciders
mycrit_ndA = [0 0 2*ones(1,myN_criteria - 2)];
mycrit_ndB = [0 0 2*ones(1,myN_criteria - 2)];

% Criteria - non-exclusive
mycrit_deciderA = [1 2 2*ones(1,myN_criteria - 2);];   % Sch A and/or B deciders
mycrit_deciderB = [2 1 2*ones(1,myN_criteria - 2);];   % Sch A and/or B deciders
mycrit_ndA = [0 2 2*ones(1,myN_criteria - 2)];
mycrit_ndB = [2 0 2*ones(1,myN_criteria - 2)];

clear mygr
i=0;

% Sch A Deciders; SchA Bndry Rel vs SchB Bndry Rel vs SchB Bndry Irrel
% i=i+1; mycrit = mycrit_multi; mygr(i) = grt; mygr(i).criteria=[mycrit]; 
% i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs=1;
% i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs=2;
% i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs=3; 
% i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs=4; 
% i=i+1; mycrit = mycrit_nondecider; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs=1; 
% i=i+1; mycrit = mycrit_nondecider; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs=2; 
% i=i+1; mycrit = mycrit_nondecider; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs=3; 
% i=i+1; mycrit = mycrit_nondecider; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs=4; 

%i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs=1; mygr(i).legend = '$Pref_{SchA}$ $Zscore(|\Delta \phi_{SchA}|)$';
i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs=1; mygr(i).legend = '$\rm \bf {Pref_{A}} \thinspace Ctg 1v2$';
i=i+1; mycrit = mycrit_ndA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs=1; mygr(i).legend = '$\rm \bf \overline{Pref_{A}} \thinspace Ctg 1v2$';
i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs=2; mygr(i).legend = '$\rm \bf {Pref_{B}} \thinspace Ctg 3v4$';
i=i+1; mycrit = mycrit_ndB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs=2; mygr(i).legend = '$\rm \bf \overline{Pref_{B}} \thinspace Ctg 3v4$';

%
% Set up sim specific parameters
sfc_mode =  22.40141110;
perm_mode = 22.40141110;

% Set stages and excludes
stage_perm = 3;
stage_sfc = 3;
opts_exclude.excludeL = 0; 
opts_exclude.excludeO = 0; 

% Load data for sp=ctg_phi and pls=ctgs
opts_pls.perm2pls_dophi=0;
opts_perm.do_phi=0;
[wrkspc_buffer, group_ctgs, ~, pls_stats_phis, bad_any_phis] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,mygr,opts_exclude,opts_pls,opts_perm);


% Load data for sp=ctg_phi and pls=ctgs
opts_pls.perm2pls_dophi=0;
opts_perm.do_phi=1;
[wrkspc_buffer, group_phi, ~, pls_stats_ctgs, bad_any_ctgs] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,mygr,opts_exclude,opts_pls,opts_perm);


%group = grouppairs_diff(group);
%
if plot_on_spect
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[group_ctgs(:)],[],opts_PM3Dcs,[]);
end

if plot_on_bargraph
    opts_PSC.hmask = blkdiag(true(2),true(2));
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group_ctgs(:),opts_PSC);
end

%
% Plot scatter
% pls_stats_ctg = -1*diff(pls_stats_ctg,[],2);
% pls_stats_ctg = pls_stats_ctg(:,1:2:end);

%
% Sch A ctgs vs Sch B Bndry
if plot_on_scatter
    %
    do_abs = 1;
    ba = bad_any_phis;
    figure;
    x = pls_stats_ctgs(~ba,1); y = pls_stats_phis(~ba,1);
    if do_abs; x=abs(x);y=abs(y); end
    [x2,y2] = reduce_xy_plusminus(x,y,opts_perm.split_plusminus);
    plot(x,y,'k.');
    hold on; plott_fit(x2,y2,'kx'); % Only fit to data points matching plusminus config.
    hold on; plot_scattergroups(x,y,[group_ctgs(1) group_phi(1)],~ba); xlabel('Ctg1-Ctg2'); ylabel('Ctg1-Ctg2 Phi');

    figure;
    x = pls_stats_ctgs(~ba,2); y = pls_stats_phis(~ba,2);
    if do_abs; x=abs(x);y=abs(y); end
    [x2,y2] = reduce_xy_plusminus(x,y,opts_perm.split_plusminus);
    plot(x,y,'k.');
    hold on; plott_fit(x2,y2,'kx'); % Only fit to data points matching plusminus config.
    hold on; plot_scattergroups(x,y,[group_ctgs(2) group_phi(2)],~ba); xlabel('Ctg3-Ctg4'); ylabel('Ctg3-Ctg4 Phi');
end


%% Fig 12 - Alpha vs Beta Scatter Sch B

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig12';

freqband1 = [10 12];
freqband2 = [16 20];

curr_sch = 'B';
switch curr_sch
    case 'A'
        ind = 1;
    case 'B'
        ind = 2;
end

% % % % % % % % % % Set up general stuff
% SFC mode
sfc_mode = 22.4014111;
perm_mode = 22.4014111;

% Stage selection
curr_stage_sp = 3;
curr_stage_sfc = 3;

opts_exclude.excludeL = 0; 
opts_exclude.excludeO = 0; 
opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC

% More pls switches
opts_pls = Opts_Pls;
opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
    opts_pls.perm2pls_dophi = 0;
    opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
    opts_pls.swap_pls = [];
    %opts_pls.swap_pls = [2,4;5,7];      
    opts_pls.do_diff = 1;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
    
    
opts_perm.do_bh0 = 1;
opts_perm.do_phi = 0;
opts_perm.split_plusminus = 0;


% % % % % % % % % % % % % % % % Pull alpha and beta pls
freqband_stats = [freqband1];
[wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
pls = out_pls.pls;
pls_stats_alpha = out_pls.pls_stats;
abscissa = out_pls.abscissa;
bad_any = out_pls.bad_any;
funames1D = out_pls.funames1D;
mypairs = out_pls.mypairs;
group0 = out_pls.group;
clear out_pls

freqband_stats = [freqband2];
[wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
pls = out_pls.pls;
pls_stats_beta = out_pls.pls_stats;
abscissa = out_pls.abscissa;
bad_any = out_pls.bad_any;
funames1D = out_pls.funames1D;
mypairs = out_pls.mypairs;
group0 = out_pls.group;
clear out_pls


% % % % % % % % % % % % % % % % Pull alpha and beta sp
freqband_stats_perm = [freqband1];
[wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm);
sp_alpha = out_perm.sig_cells;

freqband_stats_perm = [freqband2];
[wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm);
sp_beta = out_perm.sig_cells;



% % % % % % % % % % % Merge alpha and beta into single struct
sp = [sp_alpha(:,ind) sp_beta(:,ind) sp_beta(:,5)];
pls_stats = [pls_stats_alpha(:,ind) pls_stats_beta(:,ind) pls_stats_beta(:,5)];

sp = [sp_alpha(:,ind) sp_beta(:,ind) sp_beta(:,5)];
pls_stats = [pls_stats_alpha(:,ind) pls_stats_beta(:,ind) pls_stats_beta(:,5)];

%
% % % % % % % % % % % Setup group
% Create group template
mycrit = [2*ones(1,3)];
grt = Grp;
grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
clear grt2
i=0;
i=i+1; grt2(i)=grt; grt2(i).criteria(1:2)=[1 0]; grt2(i).ctgs=1; grt2(i).legend = ['Pref Sch ' curr_sch ' ' num2str(mean(freqband1)) 'Hz'];
i=i+1; grt2(i)=grt; grt2(i).criteria(1:2)=[0 1]; grt2(i).ctgs=2; grt2(i).legend = ['Pref Sch ' curr_sch ' ' num2str(mean(freqband2)) 'Hz'];
% i=i+1; grt2(i)=grt; grt2(i).criteria(5)=[1]; grt2(i).ctgs=1; grt2(i).legend = 'Rule';
i=i+1; grt2(i)=grt; grt2(i).criteria(1:2)=[1 1]; grt2(i).ctgs=1; grt2(i).legend = 'Multitasker';
group = grt2;

group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D);

%

% figure;
%     inds = 1:length(group);
%     [h1] = plot_matrix3D_custstruct([],group(inds));
    
%
 figure; % subplot(221); 
x = pls_stats(~bad_any,group(1).ctgs);
y = pls_stats(~bad_any,group(2).ctgs);
plott_fit(x,y,'k.');
hold on; plot_scattergroups(x,y,[group(:)],~bad_any);
xlabel(['FFC ' group(1).legend]); ylabel(['FFC ' group(2).legend]); % Might need fixing
xlabel(['Sch ' curr_sch ' ' num2str(mean(freqband1)) 'Hz']); ylabel(['Sch ' curr_sch ' ' num2str(mean(freqband2)) 'Hz']);



%% Figure 6E - Relevant vs irrelevant 2
% Correct figure, but doesn't work with units

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig6e';

plot_on_spect = 1;
plot_on_bargraph = 1;
plot_on_func = 0;   % Plotting within function


% SFC mode
sfc_mode = 22.4014111;
perm_mode = 22.4014111;


    % Pls switches
    freqband_stats = [16 20];
    freqband_stats_perm = [16 20];
%     freqband_stats = [10 12];
%     freqband_stats_perm = [10 12];
%     freqband_stats = [1 10];
%     freqband_stats_perm = [1 10];

    
    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 1;
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    



opts_PM3Dcs.remove_dependent = 1;
opts_PM3Dcs.stats_mode = 2;



%
% % % % Run figure for stage / animal combinations
% % % % Animal selection
% opts_exclude.excludeL = 0; 
% opts_exclude.excludeO = 1; 
% % Sample stage
% curr_stage_sp = 2; curr_stage_sfc = 2;
% [wrkspc_buffer, group_ls] = Fig_6e_rel_vs_irrel(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,opts_PM3Dcs,plot_on_func);
% % Delay selection
% curr_stage_sp = 3; curr_stage_sfc = 3;
% [wrkspc_buffer, group_ld] = Fig_6e_rel_vs_irrel(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,opts_PM3Dcs,plot_on_func);
% 
% % % % Animal selection
% opts_exclude.excludeL = 1; 
% opts_exclude.excludeO = 0; 
% % Sample stage
% curr_stage_sp = 2; curr_stage_sfc = 2;
% [wrkspc_buffer, group_os] = Fig_6e_rel_vs_irrel(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,opts_PM3Dcs,plot_on_func);
% % Delay selection
% curr_stage_sp = 3; curr_stage_sfc = 3;
% [wrkspc_buffer, group_od] = Fig_6e_rel_vs_irrel(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,opts_PM3Dcs,plot_on_func);

% % % Animal selection
opts_exclude.excludeL = 0; 
opts_exclude.excludeO = 0; 
% Sample stage
curr_stage_sp = 2; curr_stage_sfc = 2;
[wrkspc_buffer, group_bs] = Fig_6e_rel_vs_irrel(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,opts_PM3Dcs,plot_on_func);
% Delay selection
curr_stage_sp = 3; curr_stage_sfc = 3;
[wrkspc_buffer, group_bd] = Fig_6e_rel_vs_irrel(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,opts_PM3Dcs,plot_on_func);

i=0;
if plot_on_spect
    %
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[group_bs],[],opts_PM3Dcs,[]);
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[group_bd],[],opts_PM3Dcs,[]);
end

if plot_on_bargraph
    %
    opts_PSC.hmask = [];
    opts_PSC.remove_dependent = 0;
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group_bs group_bd],opts_PSC);
end

%% Figure 9B - FFC Roy Fig
% Redo Figure 9 (Roy et al Fig 4) using new generalized code
% Do this in preparation for analyzing categorization index

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig9b';


plot_on_func = 1;
% SFC mode
sfc_mode = 22.49141030;
perm_mode_sort = 22.401411100;
perm_mode_sp = 22.401411100;

% sfc_mode = 22.491510310;                 % Partial coherence all cells
% perm_mode_sort = 22.401511110;
% perm_mode_sp = 22.401511110;
% 
% sfc_mode = 22.491510310;                 % Partial coherence boundary specialist 
% perm_mode_sort = 22.401511100;
% perm_mode_sp = 22.451511100;

% sfc_mode = 41.6013111;
% perm_mode_sort = 41.6013111;
% perm_mode_sp = 41.6013111;
% sfc_mode = 52.7900000;         % Units
% sfc_mode = 52.790000000;         % Units
% perm_mode_sort = 52.700001000;   % For sorting
% perm_mode_sp = 52.700001000;     % For sp

% Setup unitsCI_vs_electBNDRY if we're using it
do_unitsCI_vs_electBNDRY = 0;

if do_unitsCI_vs_electBNDRY
    warning('FORCING: Using unit vs boundary mode!');
    sfc_mode = 52.790000001;         % Units with rel/irrel merged
    %sfc_mode = 52.790200001;         % With time series too (this is pointless for this mode)
    perm_mode_sort = 52.700001001;   % For sorting
    perm_mode_sp = 22.4514111;     % For sp 
end


% Pls switches
freqband_stats = [16 20];
freqband_stats_perm = [16 20];
% freqband_stats = [10 12];
% freqband_stats_perm = [10 12];

% Set up pls fv (popts)
% Unit exclusion
opts_exclude = Opts_Exclude;
opts_exclude.exclude_clipping = 1;
opts_exclude.exclude_60 = 0;
opts_exclude.exclude_nans = 1;
opts_exclude.excludeL = 0; 
opts_exclude.excludeO = 0; 
opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC


% More pls switches
opts_pls = Opts_Pls;
opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
    opts_pls.perm2pls_dophi = 0;
    opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
    opts_pls.swap_pls = [];
    %opts_pls.swap_pls = [2,4;5,7];      
    opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))

% Permutation test options
opts_perm = Opts_Perm;
opts_perm.do_bh0 = 1;
opts_perm.do_phi = 0;
opts_perm.split_plusminus = 0;


% Stage selection
curr_stage_sp = 2;
curr_stage_sfc = 2;

[wrkspc_buffer, mygroup, pls_sort ] = Fig_9b_Royfig_redo(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode_sort,perm_mode_sp,curr_stage_sp,freqband_stats_perm,opts_perm,do_unitsCI_vs_electBNDRY,plot_on_func);

% Stage selection
curr_stage_sp = 3;
curr_stage_sfc = 3;

[wrkspc_buffer, mygroup, pls_sort ] = Fig_9b_Royfig_redo(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode_sort,perm_mode_sp,curr_stage_sp,freqband_stats_perm,opts_perm,do_unitsCI_vs_electBNDRY,plot_on_func);


%% Figure 9C - FFC Categorization Index
% Redo Figure 9 (Roy et al Fig 4) using new generalized code
% Do this in preparation for analyzing categorization index

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig9c';


do_AB_symmetry = 0;

compare_sample_vs_delay = 1;

plot_on_func = 0;
% SFC mode
sfc_mode = 22.4914103;
perm_mode_sort = 22.4014111;
perm_mode_sp = 22.4014111;
% sfc_mode = 41.6013111;
% perm_mode_sort = 41.6013111;
% perm_mode_sp = 41.6013111;
% sfc_mode = 52.7900000;         % Units
% sfc_mode = 52.790000000;         % Units
% perm_mode_sort = 52.700001000;   % For sorting
% perm_mode_sp = 52.700001000;     % For sp



if compare_sample_vs_delay
    sfc_mode = 22.491410300;
    perm_mode_sort = 22.401411101;
    perm_mode_sp = 22.401411101;
    
    sfc_mode = 22.491510310;                 % Partial coherence all cells
    perm_mode_sort = 22.401511110;
    perm_mode_sp = 22.401511110;
    
    % sfc_mode = 41.6013111;
    % perm_mode_sort = 41.6013111;
    % perm_mode_sp = 41.6013111;
    % sfc_mode = 52.7900000;         % Units
%     sfc_mode = 52.790000001;         % Units
%     perm_mode_sort = 52.700001001;   % For sorting
%     perm_mode_sp = 52.700001001;     % For sp
end

% Setup unitsCI_vs_electBNDRY if we're using it
do_unitsCI_vs_electBNDRY = 0;

if do_unitsCI_vs_electBNDRY
    warning('FORCING: Using unit vs boundary mode!');
    sfc_mode = 52.790000001;         % Units with rel/irrel merged
    %sfc_mode = 52.790200001;         % With time series too (this is pointless for this mode)
    perm_mode_sort = 52.700001001;   % For sorting
    perm_mode_sp = 22.4514111;     % For sp    
end

% Pls switches
freqband_stats = [16 20];
freqband_stats_perm = [16 20];
% freqband_stats = [10 12];
% freqband_stats_perm = [10 12];

% Set up pls fv (popts)
% Unit exclusion
opts_exclude = Opts_Exclude;
opts_exclude.exclude_clipping = 1;
opts_exclude.exclude_60 = 0;
opts_exclude.exclude_nans = 1;
opts_exclude.excludeL = 0; 
opts_exclude.excludeO = 0; 
opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC


% More pls switches
opts_pls = Opts_Pls;
opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
    opts_pls.perm2pls_dophi = 0;
    opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
    opts_pls.swap_pls = [];
    %opts_pls.swap_pls = [2,4;5,7];      
    opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
opts_pls.collapse_pls_to_days = 0;

% Permutation test options
opts_perm = Opts_Perm;
opts_perm.do_bh0 = 1;
opts_perm.do_phi = 0;
opts_perm.split_plusminus = 0;


% Stage selection - Sample
curr_stage_sp = 2;
curr_stage_sfc = 2;

[wrkspc_buffer, mygroup, gall2s ] = Fig_9c_Cat_Index(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode_sort,perm_mode_sp,curr_stage_sp,freqband_stats_perm,opts_perm,do_unitsCI_vs_electBNDRY,plot_on_func,do_AB_symmetry);

ind = 1:length(gall2s);
% ind = [1,2];
if ~compare_sample_vs_delay
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[gall2s(ind)]);
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(gall2s(ind),opts_PSC);
end

% Stage selection - Delay
curr_stage_sp = 3;
curr_stage_sfc = 3;

[wrkspc_buffer, mygroup, gall2d ] = Fig_9c_Cat_Index(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode_sort,perm_mode_sp,curr_stage_sp,freqband_stats_perm,opts_perm,do_unitsCI_vs_electBNDRY,plot_on_func,do_AB_symmetry);


if ~compare_sample_vs_delay
    i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[gall2d(ind)]);
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(gall2d(ind),opts_PSC);
end


if compare_sample_vs_delay
    
    
    opts_PSC.do_stats = 1;
    opts_PSC.do_diagtest = 1;
    opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    opts_PSC.opts_diagtest.threshold = 0.0; 
    
    if do_AB_symmetry
        for i = 1:length(gall2s); gall2s(i).legend = ['Sample']; end
        for i = 1:length(gall2d); gall2d(i).legend = ['Delay']; end
        
        i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[gall2s(1) gall2d(1)]);    
        i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([gall2s(1) gall2d(1)],opts_PSC);
    else
        for i = 1:length(gall2s); gall2s(i).legend = [gall2s(i).legend ' Sample']; end
        for i = 1:length(gall2d); gall2d(i).legend = [gall2d(i).legend ' Delay']; end

        i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[gall2s(1) gall2d(1) gall2s(3) gall2d(3)]);
        i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([gall2s(1) gall2d(1) gall2s(3) gall2d(3)],opts_PSC);    
    end
end

%% Figure 9D - Time series analysis:
% Look for correlation between time series of unit CI and boundary response spectrogram.

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig9d';


plot_on_func = 1;
% SFC mode
sfc_mode = 23.491410101;
perm_mode_sort = 22.4014111;
perm_mode_sp = 22.4014111;
% sfc_mode = 41.6013111;
% perm_mode_sort = 41.6013111;
% perm_mode_sp = 41.6013111;
% sfc_mode = 52.7900000;         % Units
sfc_mode = 52.790200001;         % Units rel/irrel merged; timeseries
perm_mode_sort = 52.7000010;   % For sorting
perm_mode_sp = 52.7000010;     % For sp

sfc_mode2 = 23.4514111;

% Stage selection
curr_stage_sp = 2;
curr_stage_sfc = 4;

% Pls switches
freqband_stats = [16 20];
freqband_stats_perm = [16 20];
% freqband_stats = [10 12];
% freqband_stats_perm = [10 12];

% Set up pls fv (popts)
% Unit exclusion
opts_exclude = Opts_Exclude;
opts_exclude.exclude_clipping = 1;
opts_exclude.exclude_60 = 0;
opts_exclude.exclude_nans = 1;
opts_exclude.excludeL = 0; 
opts_exclude.excludeO = 0; 
opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC


% More pls switches
opts_pls = Opts_Pls;
opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
    opts_pls.perm2pls_dophi = 0;
    opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
    opts_pls.swap_pls = [];
    %opts_pls.swap_pls = [2,4;5,7];      
    opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))

opts_pls2 = opts_pls;
opts_pls2.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
opts_pls2.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    opts_pls2.do_diff = 0;
    
% Permutation test options
opts_perm = Opts_Perm;
opts_perm.do_bh0 = 1;
opts_perm.do_phi = 0;
opts_perm.split_plusminus = 0;

running_mode = 2;   % Save individual plots for boundary responses
                    % Plot correlation across all electrode pairs

[wrkspc_buffer ] = Fig_9d_time_series_plots(wrkspc_buffer,sfc_mode,sfc_mode2,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,opts_pls2,perm_mode_sort,perm_mode_sp,curr_stage_sp,freqband_stats_perm,opts_perm,running_mode);

%% Figure 9E - Time series analysis:
% Look for correlation between time series of unit CI and boundary response spectrogram.

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig9e';


plot_on_func = 1;
% SFC mode
% sfc_mode = 23.491410101;
% perm_mode_sp = 22.4014111;
% sfc_mode = 41.6013111;
% perm_mode_sort = 41.6013111;
% perm_mode_sp = 41.6013111;
% sfc_mode = 52.7900000;         % Units
sfc_mode = 52.770201001;         % Units rel/irrel merged; timeseries
perm_mode_sp = 52.770201001;     % For sp

sfc_mode2 = 23.4414111;

% Stage selection
curr_stage_sp = 2;
curr_stage_sfc = 4;

% Pls switches
freqband_stats = [16 20];
freqband_stats_perm = [16 20];
% freqband_stats = [10 12];
% freqband_stats_perm = [10 12];

% Set up pls fv (popts)
% Unit exclusion
opts_exclude = Opts_Exclude;
opts_exclude.exclude_clipping = 1;
opts_exclude.exclude_60 = 0;
opts_exclude.exclude_nans = 1;
opts_exclude.excludeL = 0; 
opts_exclude.excludeO = 1; 
opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC


% More pls switches
opts_pls = Opts_Pls;
opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
    opts_pls.perm2pls_dophi = 0;
    opts_pls.perm2pls_return_mode = 1;                 % Return mode of perm2pls (4=zscore)
    opts_pls.perm2pls_allow_signed = 0;                % 0=Abs(Diff); 1=Diff;
    opts_pls.perm2pls_split_plusminus = 2;             % 0-Return everything; 1-return pos+negative in separate columns; 2-return only positive (pls.*Ramp(Diff)); 3-return only negative (pls.*Ramp(-Diff)); 4-return only significant cells (doesn't really have a use now)
    opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
    opts_pls.swap_pls = [];
    %opts_pls.swap_pls = [2,4;5,7];      
    opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))

opts_pls2 = opts_pls;
opts_pls2.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
opts_pls2.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    opts_pls2.do_diff = 0;
    opts_pls2.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
    opts_pls2.perm2pls_allow_signed = 0;                % 0=Abs(Diff); 1=Diff;
    opts_pls2.perm2pls_split_plusminus = 0;             % 0-Return everything; 1-return pos+negative in separate columns; 2-return only positive (pls.*Ramp(Diff)); 3-return only negative (pls.*Ramp(-Diff)); 4-return only significant cells (doesn't really have a use now)
    
% Permutation test options
opts_perm = Opts_Perm;
opts_perm.do_bh0 = 1;
opts_perm.do_phi = 0;
opts_perm.split_plusminus = 0;

running_mode = 2;   % 1-Save individual plots for boundary responses
                    % 2-Plot correlation across all electrode pairs

[wrkspc_buffer ] = Fig_9e_time_series_plots_ctgsetli_mode_7(wrkspc_buffer,sfc_mode,sfc_mode2, ...
    curr_stage_sfc,freqband_stats,...
    opts_exclude,opts_pls,opts_pls2,perm_mode_sp,curr_stage_sp,freqband_stats_perm,opts_perm,running_mode);


%% Figure 13a - Does FFC boundary response enhance discriminability of FFC category response near boundary? (need to rerun with mode 22.4715111 to broaden peaks)

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig13a';


plot_on_spect = 1;
plot_on_bargraph = 1;
plot_on_func = 0;

bri_mode = 2;   % 0-Just plot traces; 1-Difference; 2-Contrast index

% SFC mode
sfc_mode =  22.4715111;
perm_mode = 22.4514111;
% sfc_mode = 41.6013111;
% % perm_mode = 41.6013111;
% sfc_mode = 52.770001001;          % For looking at unit response
% sfc_mode = 52.770201001;          % For looking at unit response

% Stage selection
curr_stage_sp = 2;
curr_stage_sfc = 2;

    % Pls switches
    freqband_stats = [16 20];
    freqband_stats_perm = [16 20];


    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 1;
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 3;  % Only take negatives
    opts_perm.do_quantiles_mode = 1;
            opts_perm.chosen_quantile = 0.1;
            opts_perm.upper_quantile = 0;
    
    
% Load pls
[wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
pls = out_pls.pls;
pls_stats = out_pls.pls_stats;
abscissa = out_pls.abscissa;
bad_any = out_pls.bad_any;
funames1D = out_pls.funames1D;
mypairs = out_pls.mypairs;
group0 = out_pls.group;
clear out_pls


% Load sp's
[wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm,opts_exclude);
sp = out_perm.sig_cells;

% Load bads perm
[bad_any_perm] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);

% Map sp's as needed
[sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any);


% Create group template
mycrit = [2*ones(1,size(sp,2))];
grt = group0(1);
grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
grt.xdata = group0(1).xdata;


% Run a simple test
clear group
i=0;
i=i+1; group(i)=grt; group(i).criteria(1)=[1]; group(i).ctgs=5;   % Ctg1-2 boundary nboundary ctgs
i=i+1; group(i)=grt; group(i).criteria(1)=[1]; group(i).ctgs=1;   % Ctg1-2 boundary noundary ctgs
i=i+1; group(i)=grt; group(i).criteria(1)=[0]; group(i).ctgs=5;   % Ctg1-2 non-deciders
i=i+1; group(i)=grt; group(i).criteria(1)=[0]; group(i).ctgs=1;   % Ctg1-2 non-deciders

i=i+1; group(i)=grt; group(i).criteria(2)=[1]; group(i).ctgs=6;   % Ctg1-2 deciders
i=i+1; group(i)=grt; group(i).criteria(2)=[1]; group(i).ctgs=2;   % Ctg1-2 deciders
i=i+1; group(i)=grt; group(i).criteria(2)=[0]; group(i).ctgs=6;   % Ctg1-2 non-deciders
i=i+1; group(i)=grt; group(i).criteria(2)=[0]; group(i).ctgs=2;   % Ctg1-2 non-deciders


% Calculate legend entries
group = group.query_legend(group0);

if bri_mode == 2
    % Calculate boundary resolvability index
    group = grouppairs_merge(group);
    pls_bnd = pls(:,:,1:4); pls_nbnd = pls(:,:,5:8);
    pls_bri = (pls_bnd - pls_nbnd) ./ (pls_nbnd + pls_bnd);
    pls = pls_bri;
    i=0;
    i=i+1; group(i).ctgs = 1; group(i).legend = 'Bnd BRI Cat vs Dog';
    i=i+1; group(i).ctgs = 1; group(i).legend = 'Non-Bnd BRI Cat vs Dog';
    i=i+1; group(i).ctgs = 2; group(i).legend = 'Bnd BRI Goc vs Tad';
    i=i+1; group(i).ctgs = 2; group(i).legend = 'Non-Bnd BRI Goc vs Tad';
end

% Load data into groups
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D);


% TAke difference
if bri_mode == 1
    group = grouppairs_merge(group);
end

% Test Plot groups, at last

if plot_on_spect
    figure;
    inds = 1:round(length(group)/2);
    [h1] = plot_matrix3D_custstruct([],group(inds));
    
    figure;
    inds = round(length(group)/2)+1:length(group);
    [h1] = plot_matrix3D_custstruct([],group(inds));
    
end

% 
opts_PSC.remove_dependent = 0;
if plot_on_bargraph
    
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
end


%% Fig 8d - Do ctg1 +ve and ctg1 -ve have different boundary responses?

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig8d';


plot_on_spect = 1;
plot_on_bargraph = 1;
plot_on_func = 0;

% SFC mode
sfc_mode = 22.4514111;
perm_mode = 22.4014111;
% sfc_mode = 41.6013111;
% perm_mode = 41.6013111;


% Stage selection
curr_stage_sp = 3;
curr_stage_sfc = 3;

    % Pls switches
    freqband_stats = [16 20];
    freqband_stats_perm = [16 20];



    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 1;
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 1;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 1;
    
        
% Load pls
[wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
pls = out_pls.pls;
pls_stats = out_pls.pls_stats;
abscissa = out_pls.abscissa;
bad_any = out_pls.bad_any;
funames1D = out_pls.funames1D;
mypairs = out_pls.mypairs;
group0 = out_pls.group;
clear out_pls

pls = pls*-1;

% Load sp's
[wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm);
sp = out_perm.sig_cells;

% Create group template
mycrit = [2*ones(1,size(sp,2))];
grt = Grp;
grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
grt.xlims_desired = [0 120]; grt.xdata = group0(1).xdata;


% Run a simple test
clear group
i=0;
i=i+1; group(i)=grt; group(i).criteria(1:4)=[1 0 0 0]; group(i).ctgs=1;
i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 0 0 0]; group(i).ctgs=1;
i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 1 0 0]; group(i).ctgs=1;
i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 0 1 0]; group(i).ctgs=2;
i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 0 0 0]; group(i).ctgs=2;
i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 0 0 1]; group(i).ctgs=2;


% Calculate legend entries
group = group.query_legend(group0);


% Load data into groups
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D);

% Test Plot groups, at last

if plot_on_spect
%     figure;
%     inds = 1:length(group);
%     [h1] = plot_matrix3D_custstruct([],group(inds));
    
    figure;
    inds = 1:3;
    [h1] = plot_matrix3D_custstruct([],group(inds));
    figure;
    inds = 4:6;
    [h1] = plot_matrix3D_custstruct([],group(inds));
end


if plot_on_bargraph
    opts_PSC.hmask = blkdiag(true(3),true(3));
    opts_PSC.remove_dependent = 0;
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
end






%% Fig 13 and Fig 13b - FFC vs Units and Units vs FFC

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);

do_units_vs_ffc = 1;    % 0-Fig13-sort by units, plot FFC (might give better sig; fewer assumptions); 1-Fig13b-sort by FFC, plot units (more intuitive)
do_bndry = 0;           % Do bndry vs non-bndry

if ~do_units_vs_ffc    
    currfigname = 'Fig13';
    
    % SFC mode
    sfc_mode = 22.4014111;
    %sfc_mode = 41.6014111;
    perm_mode = 52.7000010;
    
%     perm_mode = 22.4014111;
%     sfc_mode = 41.6014111;
else
    currfigname = 'Fig13b';

    % SFC mode
    sfc_mode = 52.7002010;
    perm_mode = 22.4014111;
    %perm_mode = 41.6014111;

end

if do_bndry
    sfc_mode = sfc_mode + 0.05
    perm_mode = perm_mode + 0.05
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
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    opts_exclude.excludeL = 0;
    opts_exclude.excludeO = 0; 
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
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
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    
    
    opts_PM3Dcs.stats_mode = 0;
    
% Plot switches
plot_on_spect = 1;
plot_on_bargraph = 1;
plot_on_func = 0;


curr_stage_sp = 2;
if do_units_vs_ffc
    curr_stage_sfc = 4;
    freqband_stats = [0.5 0.5];
else
    curr_stage_sfc = 2;
end
[wrkspc_buffer, group_bs] = Fig_13_Units_vs_Ctg(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,plot_on_func,do_units_vs_ffc)
% figure; [h1] = plot_matrix3D_custstruct([],group_bs,[],opts_PM3Dcs);


curr_stage_sp = 3;
if do_units_vs_ffc
    curr_stage_sfc = 4;
    freqband_stats = [1.5 1.5];
else
    curr_stage_sfc = 3;
end
[wrkspc_buffer, group_bd] = Fig_13_Units_vs_Ctg(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,plot_on_func,do_units_vs_ffc)

% Test Plot groups, at last
if plot_on_spect
    
    if ~do_units_vs_ffc
        figure; [h1] = plot_matrix3D_custstruct([],group_bs,[],opts_PM3Dcs);
        figure; [h1] = plot_matrix3D_custstruct([],group_bd,[],opts_PM3Dcs);
    else
        inds = 1:2; figure; [h1] = plot_matrix3D_custstruct([],group_bs(inds),[],opts_PM3Dcs);
        inds = 3:4; figure; [h1] = plot_matrix3D_custstruct([],group_bs(inds),[],opts_PM3Dcs);
        inds = 1:2; figure; [h1] = plot_matrix3D_custstruct([],group_bd(inds),[],opts_PM3Dcs);
        inds = 3:4; figure; [h1] = plot_matrix3D_custstruct([],group_bd(inds),[],opts_PM3Dcs);
%         inds = 1:2; figure; [h1] = plot_matrix3D_custstruct([],[group_bs(inds) group_bd(inds)],[],opts_PM3Dcs);
%         inds = 3:4; figure; [h1] = plot_matrix3D_custstruct([],[group_bs(inds) group_bd(inds)],[],opts_PM3Dcs);
    end
   
end

% Test bargraph

if plot_on_bargraph
    %
    
    if ~do_units_vs_ffc
        
    else
        
        sfc_mode2 = 52.7000010;
        
        if do_bndry
            sfc_mode2 = sfc_mode2 + 0.05;
        end
        
        curr_stage_sp = 2;
        curr_stage_sfc = 2;
        [wrkspc_buffer, group_bs] = Fig_13_Units_vs_Ctg(wrkspc_buffer,sfc_mode2,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,plot_on_func,do_units_vs_ffc)

        curr_stage_sp = 3;
        curr_stage_sfc = 3;
        [wrkspc_buffer, group_bd] = Fig_13_Units_vs_Ctg(wrkspc_buffer,sfc_mode2,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,plot_on_func,do_units_vs_ffc)

        
    end
    
    opts_PSC.remove_dependent=0;
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group_bs],opts_PSC);
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group_bd],opts_PSC);
end


%% Fig3c_old
% Original test code I used when determining if I would run this analysis
% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);

do_units_vs_ffc = 1;    % 0-Fig13-sort by units, plot FFC (might give better sig; fewer assumptions); 1-Fig13b-sort by FFC, plot units (more intuitive)
do_bndry = 0;           % Do bndry vs non-bndry

currfigname = 'Fig13c0';

if ~do_units_vs_ffc    
    
    
    % SFC mode
    sfc_mode = 22.4014111;
    perm_mode = 52.7000010;
else
    
    % SFC mode
    sfc_mode = 52.7002010;
    perm_mode = 22.4014111;

end

if do_bndry
    sfc_mode = sfc_mode + 0.05
    perm_mode = perm_mode + 0.05
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
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    
    
    opts_PM3Dcs.stats_mode = 0;
    
% Plot switches
plot_on_spect = 1;
plot_on_bargraph = 1;
plot_on_func = 0;


% Run analysis
opts_exclude.excludeL = 0; 
opts_exclude.excludeO = 0; 

curr_stage_sp = 2;
if do_units_vs_ffc
    curr_stage_sfc = 4;
    freqband_stats = [0.5 0.5];
else
    curr_stage_sfc = 2;
end
[wrkspc_buffer, group_bs] = Fig_13c0_Units_vs_Ctg(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,plot_on_func,do_units_vs_ffc)
% figure; [h1] = plot_matrix3D_custstruct([],group_bs,[],opts_PM3Dcs);


curr_stage_sp = 3;
if do_units_vs_ffc
    curr_stage_sfc = 4;
    freqband_stats = [1.5 1.5];
else
    curr_stage_sfc = 3;
end
[wrkspc_buffer, group_bd] = Fig_13c0_Units_vs_Ctg(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,plot_on_func,do_units_vs_ffc)

% Test Plot groups, at last
if plot_on_spect
    %
    if ~do_units_vs_ffc
        figure; [h1] = plot_matrix3D_custstruct([],group_bs,[],opts_PM3Dcs);
        figure; [h1] = plot_matrix3D_custstruct([],group_bd,[],opts_PM3Dcs);
    else
        inds = 1:3; figure; [h1] = plot_matrix3D_custstruct([],group_bs(inds),[],opts_PM3Dcs);
%         inds = 4:6; figure; [h1] = plot_matrix3D_custstruct([],group_bs(inds),[],opts_PM3Dcs);
        inds = 1:3; figure; [h1] = plot_matrix3D_custstruct([],group_bd(inds),[],opts_PM3Dcs);
%         inds = 4:6; figure; [h1] = plot_matrix3D_custstruct([],group_bd(inds),[],opts_PM3Dcs);
%         inds = 1:2; figure; [h1] = plot_matrix3D_custstruct([],[group_bs(inds) group_bd(inds)],[],opts_PM3Dcs);
%         inds = 3:4; figure; [h1] = plot_matrix3D_custstruct([],[group_bs(inds) group_bd(inds)],[],opts_PM3Dcs);
    end
   
end
%
% Test bargraph

if plot_on_bargraph
    %
    
    if ~do_units_vs_ffc
        
    else
        
        sfc_mode2 = 52.7000010;
        
        if do_bndry
            sfc_mode2 = sfc_mode2 + 0.05;
        end
        
        curr_stage_sp = 2;
        curr_stage_sfc = 2;
        [wrkspc_buffer, group_bs] = Fig_13c0_Units_vs_Ctg(wrkspc_buffer,sfc_mode2,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,plot_on_func,do_units_vs_ffc)

        curr_stage_sp = 3;
        curr_stage_sfc = 3;
        [wrkspc_buffer, group_bd] = Fig_13c0_Units_vs_Ctg(wrkspc_buffer,sfc_mode2,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,plot_on_func,do_units_vs_ffc)

        
    end
    
    opts_PSC.remove_dependent=0;
%     mask = false(6);
%     mask(1,2)=true; mask(1,3)=true; mask(4,5)=true; mask(4,6)=true;
%     opts_PSC.hmask = mask
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group_bs([1,2,1,3])],opts_PSC);
    %i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group_bs([1,3])],opts_PSC);
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group_bd([1,2,1,3])],opts_PSC);
%     i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group_bd([1,3])],opts_PSC);
end




%% Fig 13c - FFC vs Units SchA/B symmetry

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);

currfigname = 'Fig13c';

do_percent_change = 1;
do_ABsymmetry = 1;

do_units=1;             % 0-FFC; 1-Units time series; 2-units by stages

% Stage selection
curr_stage_sfc = 2;
curr_stage_sp = 2;
[wrkspc_buffer, group] = Fig_13c_Units_vs_Ctg(wrkspc_buffer,do_units,do_percent_change,do_ABsymmetry,curr_stage_sfc,curr_stage_sp);

curr_stage_sfc = 3;
curr_stage_sp = 3;
[wrkspc_buffer, group] = Fig_13c_Units_vs_Ctg(wrkspc_buffer,do_units,do_percent_change,do_ABsymmetry,curr_stage_sfc,curr_stage_sp);


do_units=2;             % 0-FFC; 1-Units time series; 2-units by stages

% Stage selection
curr_stage_sfc = 2;
curr_stage_sp = 2;
[wrkspc_buffer, group] = Fig_13c_Units_vs_Ctg(wrkspc_buffer,do_units,do_percent_change,do_ABsymmetry,curr_stage_sfc,curr_stage_sp);

curr_stage_sfc = 3;
curr_stage_sp = 3;
[wrkspc_buffer, group] = Fig_13c_Units_vs_Ctg(wrkspc_buffer,do_units,do_percent_change,do_ABsymmetry,curr_stage_sfc,curr_stage_sp);

%% Fig 13d - FFC boundary response vs FFC unit response collapse to days
% No result!!

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);

currfigname = 'Fig13d';

% Stage selection
curr_stage_sfc = 2;
curr_stage_sp = 2;

% Set up pls fv (popts)
% Unit exclusion
opts_exclude = Opts_Exclude;
opts_exclude.exclude_clipping = 1;
opts_exclude.exclude_60 = 0;
opts_exclude.exclude_nans = 1;
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 0; 
opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC



% More pls switches
opts_pls = Opts_Pls;
opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
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
opts_pls.target_pls_format = 0; % Convert pls to match this format!
opts_pls.collapse_pls_to_days = 1;

convert_to_bri = 1;

% Load Ctg response
% SFC mode
sfc_mode =  22.471411101;
perm_mode = 22.471411101;
sfc_mode  = 52.770001001;
perm_mode = 52.770001001;
[wrkspc_buffer, groupc,outc] = Fig_13d_FFC_ctgs_vs_units_bndry_regression(wrkspc_buffer, sfc_mode, perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,opts_pls,convert_to_bri)

% Load Bndry response
% SFC mode
sfc_mode =  22.451411103;
perm_mode = 22.451411103;
[wrkspc_buffer, groupb,outb] = Fig_13d_FFC_ctgs_vs_units_bndry_regression(wrkspc_buffer, sfc_mode, perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,opts_pls,convert_to_bri)

    %%
indx = 1;
indy = 1;

bad_any = outc.bad_any | outb.bad_any;
x = outc.pls_stats(~bad_any,indx);
y = outb.pls_stats(~bad_any,indy);
figure
h = plott_fit(x,y,'k.');
hold on; plot_scattergroups(x,y,[outc.group(indx), outb.group(indy)],~bad_any);


%% Fig 14b - Do units of ctg-sensitive electrode pairs have enhanced irrelevant attenuation (no real differences)

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig14b';

plot_on_spect = 1;
plot_on_bargraph = 1;

% SFC mode
% sfc_mode = 22.4014111;
sfc_mode = 52.7003010;
perm_mode = 22.4014111
% perm_mode = 52.7000010;

% Stage selection
curr_stage_sfc = 4;
curr_stage_sp = 2;


    % Pls switches
    freqband_stats = [1.6, 1.6];
    freqband_stats_perm = [16 20];


    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 1;
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    
    warning('To do - compare sch sensitive electrodes ability to attenuate irrelevant ctg response to other electrodes.');

        
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

[wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm);
sp = out_perm.sig_cells;
[sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md);

% Create group template
mycrit = [2*ones(1,size(sp,2))];
grt = group0(1);
grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;



clear group
i=0;
if opts_perm.split_plusminus == 0 || opts_perm.split_plusminus == 2 
    i=i+1; group(i)=grt; group(i).criteria(1)=[1]; group(i).ctgs=1;   % Ctg1-2 deciders
    i=i+1; group(i)=grt; group(i).criteria(1)=[1]; group(i).ctgs=3;   % Ctg1-2 non-deciders
    i=i+1; group(i)=grt; group(i).criteria(1)=[0]; group(i).ctgs=1;   % Ctg1-2 deciders
    i=i+1; group(i)=grt; group(i).criteria(1)=[0]; group(i).ctgs=3;   % Ctg1-2 non-deciders
end
if opts_perm.split_plusminus == 0 || opts_perm.split_plusminus == 3 
    i=i+1; group(i)=grt; group(i).criteria(2)=[1]; group(i).ctgs=2;   % Ctg1-2 deciders
    i=i+1; group(i)=grt; group(i).criteria(2)=[1]; group(i).ctgs=4;   % Ctg1-2 non-deciders            i=i+1; group(i)=grt; group(i).criteria(5)=[0]; group(i).ctgs=2;   % Ctg1-2 deciders
    i=i+1; group(i)=grt; group(i).criteria(2)=[0]; group(i).ctgs=2;   % Ctg1-2 non-deciders
    i=i+1; group(i)=grt; group(i).criteria(2)=[0]; group(i).ctgs=4;   % Ctg1-2 non-deciders
end

% Calculate legend entries
group = group.query_legend(group0);


% Load data into groups
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);

%group = grouppairs_diff(group);

if plot_on_spect
    
    if length(group) == 8
        figure;
        inds = 1:4;
        [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D);
        figure;
        inds = 5:8;
        [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D);
    else
        figure;
        inds = 1:length(group);
        [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D);
    end
end

if plot_on_bargraph
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PM3Dcs);
end






%% Fig 8b redo (Fig 8b2) - Ctgmode 0 vs mode 5 with scatterplot, units
% % Found no clear differences between FFC and units in this case (both ffc
% % and units multitask ctgs and boundaries to some extent)
% Set up default parameters
% Has 2 modes, defined by Fig8b2_mode
% 1) Plot spectrum
% 2) Plot scatter
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig8b2';

Fig8b2_mode = 1;

plot_on_spect = 0;
plot_on_scatter = 0;
switch Fig8b2_mode
    case 1; plot_on_spect=1;
    case 2; plot_on_scatter = 1;
end

do_units=1;

if ~do_units
    % SFC mode FFC
    sfc_mode = 22.4014111;
    perm_mode = 22.4014111;
    sfc_mode2 = 22.4514111;
    perm_mode2 = 22.4514111;
    pername = 'FFC';
else
    % SFC mode units
    sfc_mode = 52.7003010;
    perm_mode = 52.7000010;
    sfc_mode2 = 52.7502010;
    perm_mode2 = 52.7500010;
    pername = 'Units';
end

% Stage selection
curr_stage_sfc = 4;
curr_stage_sp = 3;


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
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    
    % Plotting options
    paperfig_mode = 0;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    
    % Load pls
    [wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
    pls1 = out_pls.pls;
    pls_stats1 = out_pls.pls_stats;
    abscissa = out_pls.abscissa;
    abscissa2 = out_pls.abscissa2;
    bad_any1 = out_pls.bad_any;
    funames1D = out_pls.funames1D;
    mypairs = out_pls.mypairs;
    group01 = out_pls.group;
    clear out_pls
    
    [wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode2,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
    pls2 = out_pls.pls;
    pls_stats2 = out_pls.pls_stats;
    bad_any2 = out_pls.bad_any;
    group02 = out_pls.group;
    clear out_pls
    
    pls = cat(3,pls1,pls2);
    pls_stats = cat(2,pls_stats1,pls_stats2);
    bad_any = bad_any1 | bad_any2;
    group0 = [group01 group02];
    
    for i=length(group01)+1:length(group0)          % Fix the ctg indicies to correspond to pls values; this is important for getting legend entries later
        group0(i).ctgs = group0(i).ctgs + length(group01);  
    end
    
    clear pls1 pls2 pls_stats1 pls_stats2 bad_any1 bad_any2 group01 group02
%
    % Group based on perm_mode
    % Load sp's
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm);
    sp1 = out_perm.sig_cells;
    
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode2,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm);
    sp2 = out_perm.sig_cells;
    
    sp = cat(2,sp1,sp2);
    clear sp1 sp2

    % Load bads perm
    [bad_any_perm1] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);
    [bad_any_perm2] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);
    bad_any_perm = bad_any_perm1 | bad_any_perm2;
    clear bad_any_perm1 bad_any_perm2

    % Map sp's as needed
    [sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any);

    
    %
% Create group template
mycrit = [2*ones(1,size(sp,2))];
grt = group0(1);
grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;

switch Fig8b2_mode
    case 1 
        % Spectras
        clear group
        i=0;
        i=i+1; group(i)=grt; group(i).criteria([1,6])=[1,2]; group(i).ctgs=6; group(i).legend = ['$\rm \bf Ctgs \thinspace C/D \thinspace ' pername '  $']; group(i).data_name = group0(6).legend;
        i=i+1; group(i)=grt; group(i).criteria([1,6])=[0,2]; group(i).ctgs=6; group(i).legend = ['$\rm \bf Ctgs \thinspace \overline{C/D} \thinspace ' pername '  $']; group(i).data_name = group0(6).legend;
        i=i+1; group(i)=grt; group(i).criteria([2,7])=[1,2]; group(i).ctgs=7; group(i).legend = ['$\rm \bf Ctgs \thinspace F/T \thinspace ' pername '  $'];  group(i).data_name = group0(7).legend;
        i=i+1; group(i)=grt; group(i).criteria([2,7])=[0,2]; group(i).ctgs=7; group(i).legend = ['$\rm \bf Ctgs \thinspace \overline{F/T} \thinspace ' pername '  $'];  group(i).data_name = group0(7).legend;
        
        for i = 1:4; group(i).interpreter = 'latex'; end
        % Calculate legend entries
        %group = group.query_legend(group0);
    case 2
        % Scatterplots
        clear group
        i=0;
        i=i+1; group(i)=grt; group(i).criteria([1,6])=[1,0]; group(i).ctgs=1;   
        i=i+1; group(i)=grt; group(i).criteria([1,6])=[0,1]; group(i).ctgs=6;   
        i=i+1; group(i)=grt; group(i).criteria([1,6])=[1,1]; group(i).ctgs=1;   
        i=i+1; group(i)=grt; group(i).criteria([2,7])=[1,0]; group(i).ctgs=2;   
        i=i+1; group(i)=grt; group(i).criteria([2,7])=[0,1]; group(i).ctgs=7;   
        i=i+1; group(i)=grt; group(i).criteria([2,7])=[1,1]; group(i).ctgs=2;
        % Calculate legend entries
        group = group.query_legend(group0);
        group(3).legend = 'Ctgs & Bndry';
        group(6).legend = 'Ctgs & Bndry';
end



% Load data into groups
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);

%
% Test Plot groups, at last
if plot_on_spect
    opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
    
    inds = 1:2;
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    
    inds = 3:4;
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    
end

if plot_on_scatter
    % plot_on_scatter
    figure; % subplot(221); 
    ind1 = 1;
    ind2 = 2;
    ind3 = 3;
    x = pls_stats(~bad_any,group(ind1).ctgs);
    y = pls_stats(~bad_any,group(ind2).ctgs);
    plott_fit(x,y,'k.');
    hold on; plot_scattergroups(x,y,[group(ind1),group(ind2),group(ind3)],~bad_any);
    xlabel(['FFC ' group(ind1).legend]); ylabel(['FFC ' group(ind2).legend]); % Might need fixing
    
    
    fprintf(['L1 Slope SchA Rel Vs Irrel= ' num2str(std(pls_stats(~bad_any,3))/std(pls_stats(~bad_any,1))) '\n']);
    fprintf(['L1 Slope SchB Rel Vs Irrel= ' num2str(std(pls_stats(~bad_any,4))/std(pls_stats(~bad_any,2))) '\n']);
    
    figure; % subplot(221); 
    ind1 = 4;
    ind2 = 5;
    ind3 = 6;
    x = pls_stats(~bad_any,group(ind1).ctgs);
    y = pls_stats(~bad_any,group(ind2).ctgs);
    plott_fit(x,y,'k.');
    hold on; plot_scattergroups(x,y,[group(ind1),group(ind2),group(ind3)],~bad_any);
    xlabel(['FFC ' group(ind1).legend]); ylabel(['FFC ' group(ind2).legend]); % Might need fixing
    
    
    fprintf(['L1 Slope SchA Rel Vs Irrel= ' num2str(std(pls_stats(~bad_any,3))/std(pls_stats(~bad_any,1))) '\n']);
    fprintf(['L1 Slope SchB Rel Vs Irrel= ' num2str(std(pls_stats(~bad_any,4))/std(pls_stats(~bad_any,2))) '\n']);
    
end




%% Fig 8c redo (Fig 8c2) - Ctgmode 0 vs mode 5 with units
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig8c2';

plot_on_spect = 1;
plot_on_bargraph = 1;

do_units=0;

pername='';
if ~do_units
    % SFC mode FFC
    sfc_mode = 22.4014111;
    perm_mode = 22.4014111;
    sfc_mode2 = 22.4514111;
    perm_mode2 = 22.4514111;
    %pername = 'FFC';
else
    % SFC mode units
    sfc_mode = 52.7003010;
    perm_mode = 52.7000010;
    sfc_mode2 = 52.7503010;
    perm_mode2 = 52.7500010;
    %pername = 'Units';
end

% Stage selection
curr_stage_sfc = 2;
curr_stage_sp = 2;

if do_units
    curr_stage_sfc=4;
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
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    
    % Plotting options
    paperfig_mode = 0;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    
    % Load pls
    [wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
    pls1 = out_pls.pls;
    pls_stats1 = out_pls.pls_stats;
    abscissa = out_pls.abscissa;
    abscissa2 = out_pls.abscissa2;
    bad_any1 = out_pls.bad_any;
    funames1D = out_pls.funames1D;
    mypairs = out_pls.mypairs;
    group01 = out_pls.group;
    clear out_pls
    
    [wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode2,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
    pls2 = out_pls.pls;
    pls_stats2 = out_pls.pls_stats;
    bad_any2 = out_pls.bad_any;
    group02 = out_pls.group;
    clear out_pls
    
    pls = cat(3,pls1,pls2);
    pls_stats = cat(2,pls_stats1,pls_stats2);
    bad_any = bad_any1 | bad_any2;
    group0 = [group01 group02];
    
    
    if opts_pls.do_diff || opts_pls.perm2pls_allow_signed     % If taking difference, reverse for boundary
        pls=pls*-1; pls_stats=pls_stats*-1;
    end
    
    for i=length(group01)+1:length(group0)          % Fix the ctg indicies to correspond to pls values; this is important for getting legend entries later
        group0(i).ctgs = group0(i).ctgs + length(group01);  
    end
    
    clear pls1 pls2 pls_stats1 pls_stats2 bad_any1 bad_any2 group01 group02
%
    % Group based on perm_mode
    % Load sp's
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm);
    sp1 = out_perm.sig_cells;
    
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode2,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm);
    sp2 = out_perm.sig_cells;
    
    sp = cat(2,sp1,sp2);
    clear sp1 sp2

    % Load bads perm
    [bad_any_perm1] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);
    [bad_any_perm2] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);
    bad_any_perm = bad_any_perm1 | bad_any_perm2;
    clear bad_any_perm1 bad_any_perm2

    % Map sp's as needed
    [sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any);

    
    %
% Create group template
mycrit = [2*ones(1,size(sp,2))];
grt = group0(1);
grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;


% Spectras
clear group
i=0;
i=i+1; group(i)=grt; group(i).criteria([1,2,6,7])=[1,2,2,2]; group(i).ctgs=6; group(i).legend = ['Cat/Dog Bndry' pername ''];
i=i+1; group(i)=grt; group(i).criteria([1,2,6,7])=[1,2,2,2]; group(i).ctgs=7; group(i).legend = ['Goc/Tad Bndry' pername ''];
i=i+1; group(i)=grt; group(i).criteria([1,2,6,7])=[2,1,2,2]; group(i).ctgs=7; group(i).legend = ['Goc/Tad Bndry' pername ''];
i=i+1; group(i)=grt; group(i).criteria([1,2,6,7])=[2,1,2,2]; group(i).ctgs=6; group(i).legend = ['Cat/Dog Bndry' pername ''];

% % Spectras
% clear group
% i=0;
% i=i+1; group(i)=grt; group(i).criteria([1,2,6,7])=[1,0,2,2]; group(i).ctgs=6; group(i).legend = ['Cat/Dog Pref' pername ''];
% i=i+1; group(i)=grt; group(i).criteria([1,2,6,7])=[0,1,2,2]; group(i).ctgs=6; group(i).legend = ['Goc/Tad Pref' pername ''];
% i=i+1; group(i)=grt; group(i).criteria([1,2,6,7])=[0,0,2,2]; group(i).ctgs=6; group(i).legend = ['Neither' pername ''];
% i=i+1; group(i)=grt; group(i).criteria([1,2,6,7])=[1,0,2,2]; group(i).ctgs=7; group(i).legend = ['Cat/Dog Pref' pername ''];
% i=i+1; group(i)=grt; group(i).criteria([1,2,6,7])=[0,1,2,2]; group(i).ctgs=7; group(i).legend = ['Goc/Tad Pref' pername ''];
% i=i+1; group(i)=grt; group(i).criteria([1,2,6,7])=[0,0,2,2]; group(i).ctgs=7; group(i).legend = ['Neither' pername ''];


for i = 1:length(group); group(i).data_name = '\Delta FFC Bnd NBnd'; end
%for i = 1:4; group(i).interpreter = 'latex'; end
% Calculate legend entries
%group = group.query_legend(group0);

% Load data into groups
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);

%
% Test Plot groups, at last
if plot_on_spect
    opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
    
    inds = 1:2;
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    
    inds = 3:4;
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    
end

if plot_on_bargraph
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group);
end

%% Fig 8c redo2 (Fig 8c3) - Ctgmode 0 vs mode 5 with units - allow for % changes
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig8c3';

plot_on_spect = 1;
plot_on_bargraph = 1;

do_units=0;             % 0-FFC; 1-Units time series; 2-units by stages
do_percent_change = 1;

pername='';
if do_units == 0            % FFC
    % SFC mode FFC
    sfc_mode = 22.451411103;
    perm_mode = 22.401411101;
    %pername = 'FFC';
elseif do_units == 1        % Units time series
    % SFC mode units
    sfc_mode = 52.7502010;
    perm_mode = 52.7000010;
    %pername = 'Units';
elseif do_units == 2        % Units by stages
    sfc_mode = 52.7500010;
    perm_mode = 52.7000010;    
end

% Stage selection
curr_stage_sfc = 2;
curr_stage_sp = 2;

if do_units == 1
    curr_stage_sfc=4;
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
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 1;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;   % Note when I first submit the paper, this was set to 0 for some reason!
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    opts_perm.alpha_bh0 = 0.2;
%     opts_perm.alpha_bh0 = 0.05;
    
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dcs.remove_dependent=0;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    
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
    
    
    if opts_pls.do_diff || opts_pls.perm2pls_allow_signed     % If taking difference, reverse for boundary
        pls=pls*-1; pls_stats=pls_stats*-1;
    end
    
    
    % pls for % change (if turned on)
    if do_percent_change
        opts_pls.permdat2pls = 1;
        opts_pls.perm2pls = 0;
        opts_pls.do_diff = 0;
        [wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
        pls2 = out_pls.pls;
        pls_stats2 = out_pls.pls_stats;
        % Make same dimensionality as other pls
        pls2 = pls2(:,:,1:2:end); 
        pls_stats2 = pls_stats2(:,1:2:end); 
        clear out_pls
    end
    
    % Group based on perm_mode
    % Load sp's
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm,opts_exclude);
    sp = out_perm.sig_cells;
    
    % Load bads perm
    [bad_any_perm] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);
 
    % Map sp's as needed
    [sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any);


% Create group template
mycrit = [2*ones(1,size(sp,2))];
grt = group0(1);
grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
%grt.ylims_desired = [-20 20];

% Spectras
clear group
i=0;
i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,0]; group(i).ctgs=1;
i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,0]; group(i).ctgs=2;
i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,1]; group(i).ctgs=1;
i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,1]; group(i).ctgs=2;

% % Spectras
% clear group
% i=0;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,0]; group(i).ctgs=1;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,0]; group(i).ctgs=1;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,1]; group(i).ctgs=2;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,0]; group(i).ctgs=2;

% % Spectras
% clear group
% i=0;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,0]; group(i).ctgs=1;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,0]; group(i).ctgs=1;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,1]; group(i).ctgs=1;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,0]; group(i).ctgs=2;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,0]; group(i).ctgs=2;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,1]; group(i).ctgs=2;

% Load data into groups
if do_percent_change
    group2 = get_grouped_cells_and_data(group,sp,pls2,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);
end
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);


if do_percent_change
    group_temp = [group; group2];   % Alternate group and group2 entries
    group_temp = group_temp(:)';
    group = grouppairs_merge(group_temp,opts_pls.perm2pls_dophi,2);
    %for i = 1:length(group);group(i).data_name = '% change'; end
end

i=0;
i=i+1; group(i).legend = 'Cat/Dog Bndry';
% i=i+1; group(i).legend = 'No pref';
i=i+1; group(i).legend = 'Goc/Tad Bndry';
i=i+1; group(i).legend = 'Cat/Dog Bndry';
i=i+1; group(i).legend = 'Goc/Tad Bndry';
% i=i+1; group(i).legend = 'No pref';




% Test Plot groups, at last
if plot_on_spect
    
    opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
    
    inds = 1:2
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    
    inds = 3:4
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);

    
end

if plot_on_bargraph
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
    if do_units >= 1; title = ''; end
end


%% Fig 8c redo2 (Fig 8c4) - Ctgmode 0 vs mode 5 - compare electrodes instead of ctgs
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig8c4';

plot_on_spect = 1;
plot_on_bargraph = 1;

do_units=0;             % 0-FFC; 1-Units time series; 2-units by stages
do_percent_change = 1;

group_pref_non_pref = 0;

pername='';
if do_units == 0            % FFC
    % SFC mode FFC
    sfc_mode = 22.451411103;
    perm_mode = 22.401411101;
    %pername = 'FFC';
elseif do_units == 1        % Units time series
    % SFC mode units
    sfc_mode = 52.7502010;
    perm_mode = 52.7000010;
    %pername = 'Units';
elseif do_units == 2        % Units by stages
    sfc_mode = 52.7500010;
    perm_mode = 52.7000010;    
end

% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 3;

if do_units == 1
    curr_stage_sfc=4;
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
    opts_exclude.exclude_nans = 0;
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 1;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;   % Note when I first submit the paper, this was set to 0 for some reason!
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    opts_perm.alpha_bh0 = 0.2;
%     opts_perm.alpha_bh0 = 0.05;
    
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dcs.remove_dependent = 0;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    
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
    
    
    if opts_pls.do_diff || opts_pls.perm2pls_allow_signed     % If taking difference, reverse for boundary
        pls=pls*-1; pls_stats=pls_stats*-1;
    end
    
    
    % pls for % change (if turned on)
    if do_percent_change
        opts_pls.permdat2pls = 1;
        opts_pls.perm2pls = 0;
        opts_pls.do_diff = 0;
        [wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
        pls2 = out_pls.pls;
        pls_stats2 = out_pls.pls_stats;
        % Make same dimensionality as other pls
        pls2 = pls2(:,:,1:2:end); 
        pls_stats2 = pls_stats2(:,1:2:end); 
        clear out_pls
    end
    
    % Group based on perm_mode
    % Load sp's
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm,opts_exclude);
    sp = out_perm.sig_cells;
    
    % Load bads perm
    [bad_any_perm] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);
 
    % Map sp's as needed
    [sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any);


% Create group template
mycrit = [2*ones(1,size(sp,2))];
grt = group0(1);
grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
%grt.ylims_desired = [-20 20];

% Spectras
if ~group_pref_non_pref
    clear group
    i=0;
    i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,0]; group(i).ctgs=1;
    i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,1]; group(i).ctgs=1;
    i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,0]; group(i).ctgs=1;
    i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,0]; group(i).ctgs=2;
    i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,1]; group(i).ctgs=2;
    i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,0]; group(i).ctgs=2;

else

    clear group
    i=0;
    i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,1]; group(i).ctgs=1;
    i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,0]; group(i).ctgs=1;
    i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,1]; group(i).ctgs=1;
    i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,0]; group(i).ctgs=1;
    i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,1]; group(i).ctgs=2;
    i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,0]; group(i).ctgs=2;
    i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,1]; group(i).ctgs=2;
    i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,0]; group(i).ctgs=2;

    group2 = [group_merge(group(1:3)), group(4), group_merge(group(5:7)), group(8)];
    group = group2; clear group2

end

% Load data into groups
if do_percent_change
    group2 = get_grouped_cells_and_data(group,sp,pls2,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);
end
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);


if do_percent_change
    group_temp = [group; group2];   % Alternate group and group2 entries
    group_temp = group_temp(:)';
    group = grouppairs_merge(group_temp,opts_pls.perm2pls_dophi,2);
    %for i = 1:length(group);group(i).data_name = '% change'; end
end


if ~group_pref_non_pref
    i=0;
    i=i+1; group(i).legend = 'Cat/Dog Pref';
    i=i+1; group(i).legend = 'Goc/Tad Pref';
    i=i+1; group(i).legend = 'No pref';
    i=i+1; group(i).legend = 'Cat/Dog Pref';
    i=i+1; group(i).legend = 'Goc/Tad Pref';
    i=i+1; group(i).legend = 'No pref';
else
    i=0;
    i=i+1; group(i).legend = 'CS';
    i=i+1; group(i).legend = 'Not CS';
    i=i+1; group(i).legend = 'CS';
    i=i+1; group(i).legend = 'Not CS';
end



% Test Plot groups, at last
if plot_on_spect
    
    opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
    
    inds = 1:3;
    if group_pref_non_pref
        inds = 1:2; end
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    
    inds = 4:6;
    if group_pref_non_pref
        inds = 3:4; end
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);

    
end

if plot_on_bargraph
    group = group.arrayload('zlims_desired',[0 18]);
    if ~group_pref_non_pref
%         opts_PSC.hmask = blkdiag(ones(3),ones(3));
%         i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
%         if do_units >= 1; title = ''; end
    
    
        opts_PSC.hmask = blkdiag(ones(3));
        %opts_PSC.hmask = blkdiag(zeros(3)); opts_PSC.hmask(:,end) = 1;
        i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group(1:3),opts_PSC);
        i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group(4:6),opts_PSC);
    else
        i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group(:),opts_PSC);
    end
end

%%
%% Fig 8e - (Based on Fig8c2) - Ctgmode 7 vs mode 5 with units - allow for % changes
% Instead do 60% and 100% deciders individually
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig8e';

plot_on_spect = 1;
plot_on_bargraph = 1;

do_units=0;             % 0-FFC; 1-Units time series; 2-units by stages
do_percent_change = 1;

pername='';
if do_units == 0            % FFC
    % SFC mode FFC
    sfc_mode = 22.451411103;
    perm_mode = 22.471411101;
    %pername = 'FFC';
elseif do_units == 1        % Units time series
    % SFC mode units
    sfc_mode = 52.7502010;
    perm_mode = 52.7000010;
    %pername = 'Units';
elseif do_units == 2        % Units by stages
    sfc_mode = 52.7500010;
    perm_mode = 52.7000010;    
end

% Stage selection
curr_stage_sfc = 2;
curr_stage_sp = 2;

if do_units == 1
    curr_stage_sfc=4;
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
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 1;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 0;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    
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
    
    
    if opts_pls.do_diff || opts_pls.perm2pls_allow_signed     % If taking difference, reverse for boundary
        pls=pls*-1; pls_stats=pls_stats*-1;
    end
    
    
    % pls for % change (if turned on)
    if do_percent_change
        opts_pls.permdat2pls = 1;
        opts_pls.perm2pls = 0;
        opts_pls.do_diff = 0;
        [wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
        pls2 = out_pls.pls;
        pls_stats2 = out_pls.pls_stats;
        % Make same dimensionality as other pls
        pls2 = pls2(:,:,1:2:end); 
        pls_stats2 = pls_stats2(:,1:2:end); 
        clear out_pls
    end
    
    % Group based on perm_mode
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


% Spectras
clear group
ctg0 = 1;
i=0;
i=i+1; group(i)=grt; group(i).criteria([1,2,5,6])=[1,1,0,0]; group(i).ctgs=ctg0;
i=i+1; group(i)=grt; group(i).criteria([1,2,5,6])=[1,0,0,0]; group(i).ctgs=ctg0;
i=i+1; group(i)=grt; group(i).criteria([1,2,5,6])=[0,1,0,0]; group(i).ctgs=ctg0;

i=i+1; group(i)=grt; group(i).criteria([1,2,5,6])=[0,0,1,1]; group(i).ctgs=ctg0;
i=i+1; group(i)=grt; group(i).criteria([1,2,5,6])=[0,0,1,0]; group(i).ctgs=ctg0;
i=i+1; group(i)=grt; group(i).criteria([1,2,5,6])=[0,0,0,1]; group(i).ctgs=ctg0;

i=i+1; group(i)=grt; group(i).criteria([1,2,5,6])=[0,0,0,0]; group(i).ctgs=ctg0;

group3 = group_merge(group(1:3));
group4 = group_merge(group(4:6));

group = [group3, group4, group(end)];



% i=i+1; group(i)=grt; group(i).criteria([1,5])=[1,0]; group(i).ctgs=2;
% i=i+1; group(i)=grt; group(i).criteria([1,5])=[0,1]; group(i).ctgs=2;
% i=i+1; group(i)=grt; group(i).criteria([2,6])=[1,0]; group(i).ctgs=2;
% i=i+1; group(i)=grt; group(i).criteria([2,6])=[0,1]; group(i).ctgs=2;


% % Spectras
% clear group
% i=0;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,0]; group(i).ctgs=1;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,0]; group(i).ctgs=1;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,1]; group(i).ctgs=2;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,0]; group(i).ctgs=2;

% % Spectras
% clear group
% i=0;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,0]; group(i).ctgs=1;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,0]; group(i).ctgs=1;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,1]; group(i).ctgs=1;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,0]; group(i).ctgs=2;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,0]; group(i).ctgs=2;
% i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,1]; group(i).ctgs=2;

% Load data into groups
if do_percent_change
    group2 = get_grouped_cells_and_data(group,sp,pls2,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);
end
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);


if do_percent_change
    group_temp = [group; group2];   % Alternate group and group2 entries
    group_temp = group_temp(:)';
    group = grouppairs_merge(group_temp,opts_pls.perm2pls_dophi,2);
    %for i = 1:length(group);group(i).data_name = '% change'; end
end


i=0;
i=i+1; group(i).legend = '60% deciders';
i=i+1; group(i).legend = '100% deciders';
i=i+1; group(i).legend = 'non-deciders';

% i=0;
% i=i+1; group(i).legend = 'Cat/Dog Bndry';
% % i=i+1; group(i).legend = 'No pref';
% i=i+1; group(i).legend = 'Goc/Tad Bndry';
% i=i+1; group(i).legend = 'Goc/Tad Bndry';
% i=i+1; group(i).legend = 'Cat/Dog Bndry';
% % i=i+1; group(i).legend = 'No pref';




% Test Plot groups, at last
if plot_on_spect
    %%
    opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
    
    inds = 1:length(group);
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    
%     inds = 5:8
%     figure;
%     [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);

    
end

if plot_on_bargraph
    %%
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
    if do_units >= 1; title = ''; end
end



%% Fig 8a redo (Fig 8a2) - Ctgmode 0 vs mode 4 with units
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig8a2';

plot_on_spect = 1;
plot_on_bargraph = 1;

do_units=0;             % 0-FFC; 1-Units time series; 2-units by stages
do_difference = 1;
do_percent_change = 1;

pername='';
if do_units == 0            % FFC
    % SFC mode FFC
    sfc_mode = 22.4414111;
    perm_mode = 22.4014111;
    %pername = 'FFC';
elseif do_units == 1        % Units time series
    % SFC mode units
    sfc_mode = 52.7402010;
    perm_mode = 52.7000010;
    %pername = 'Units';
elseif do_units == 2        % Units by stages
    sfc_mode = 52.7400010;
    perm_mode = 52.7000010;    
end

% Stage selection
curr_stage_sfc = 2;
curr_stage_sp = 2;

if do_units == 1
    curr_stage_sfc=4;
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
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    
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
  
    
    if opts_pls.do_diff || opts_pls.perm2pls_allow_signed     % If taking difference, reverse for boundary
        pls=pls*-1; pls_stats=pls_stats*-1;
    end

    
    clear pls1 pls2 pls_stats1 pls_stats2 bad_any1 bad_any2 group01 group02
%
    % Group based on perm_mode
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

% Spectras
clear group
i=0;
i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,0]; group(i).ctgs=6;
i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,1]; group(i).ctgs=6;
i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,1]; group(i).ctgs=6;
i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,0]; group(i).ctgs=5;
i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,1]; group(i).ctgs=5;
i=i+1; group(i)=grt; group(i).criteria([1,2])=[1,1]; group(i).ctgs=5;
i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,0]; group(i).ctgs=6;
i=i+1; group(i)=grt; group(i).criteria([1,2])=[0,0]; group(i).ctgs=5;

group = group.query_legend(group0);

group2(1) = group_merge(group(1:3));
group2(2) = group_merge(group(4:6));
group2(3) = group(7);
group2(4) = group(8);
group=group2;


% Load data into groups
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);

if do_difference
    group = grouppairs_merge(group,opts_pls.perm2pls_dophi,do_percent_change);
    
    for i = 1:length(group); group(i).data_name = '% change'; end

    i=0;
    i=i+1; group(i).legend = 'Ens.';
    i=i+1; group(i).legend = 'Non Ens.';
end

% Test Plot groups, at last
if plot_on_spect
    opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
    
    inds = 1:length(group);
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
%     
%     inds = 3:4;
%     figure;
%     [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    
end

if plot_on_bargraph
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
    if do_units >= 1; title = ''; end
end


%% Fig 4C2 - Redo Boundary Ctgsetli mode 4
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig4c2';

plot_on_spect = 1;
plot_on_bargraph = 1;
plot_on_func = 0;

do_units = 2;

% SFC mode
if do_units ==0 
    sfc_mode = 22.441411102;
    perm_mode = 22.441411102;
    
%     sfc_mode =  22.451511113; % Partial
%     sfc_mode =  22.441511112; % Partial

    % sfc_mode = 2.3015111;
    % perm_mode = 2.3015111;
    % perm_mode = 5;
    % sfc_mode = 41.6514111;
    % perm_mode = 41.6514111;
elseif do_units == 1
    sfc_mode = 52.740201002;
    perm_mode = 52.740001002;
elseif do_units == 2
    sfc_mode = 52.740001002;
    perm_mode = 52.740001002;
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
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    
    [wrkspc_buffer, group, groupm ] = Fig_4c2_Bndry_redo(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,opts_PM3Dcs,opts_PM3Dsp,plot_on_func)

    inds1=[6,5,8,7,2,1,4,3];
    inds2=[3,4,1,2];
    
    if plot_on_spect
        opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
        group_pl = group;
        
        shiftinds = 3:4
        shiftval = -0.1;
        inds = [6,5,8,7];
        group_pl(inds(shiftinds)) = group_shift_data(group_pl(inds(shiftinds)), shiftval);
        figure;
        [h1] = plot_matrix3D_custstruct([],group_pl(inds),opts_PM3D,opts_PM3Dcs);
        
        
        shiftinds = 3:4
        shiftval = -0.1;
        inds = [2,1,4,3];
        group_pl(inds(shiftinds)) = group_shift_data(group_pl(inds(shiftinds)), shiftval);
        figure;
        [h1] = plot_matrix3D_custstruct([],group_pl(inds),opts_PM3D,opts_PM3Dcs);
        
        figure; 
        [h1] = plot_matrix3D_custstruct([],groupm(inds2),opts_PM3D,opts_PM3Dcs);

    end
    

    % Test bargraph
    if plot_on_bargraph
        fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group(inds1),opts_PSC);
        fign; [h1, h, p, out.PSC] = plot_stats_custstruct(groupm(inds2),opts_PSC);
    end
    
%% Fig 4C3 - Redo Boundary Ctgsetli mode 4
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig4c3';

plot_on_spect = 1;
plot_on_bargraph = 1;
plot_on_func = 0;

do_units = 0;

% SFC mode
if do_units ==0 
    sfc_mode = 22.441311102;
    perm_mode = 22.441311102;
    
%     sfc_mode =  22.451511113; % Partial
%     sfc_mode =  22.441511112; % Partial

    % sfc_mode = 2.3015111;
    % perm_mode = 2.3015111;
    % perm_mode = 5;
    % sfc_mode = 41.6514111;
    % perm_mode = 41.6514111;
elseif do_units == 1
    sfc_mode = 52.740201002;
    perm_mode = 52.740001002;
elseif do_units == 2
    sfc_mode = 52.740001002;
    perm_mode = 52.740001002;
end



if do_units == 1; curr_stage_sfc = 4; end

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
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    
    
    % Stage selection
    curr_stage_sfc = 2;
    curr_stage_sp = 2;
    if do_units == 1; curr_stage_sfc = 4; end
    [wrkspc_buffer, group_s, groupm_s ] = Fig_4c2_Bndry_redo(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,opts_PM3Dcs,opts_PM3Dsp,plot_on_func)
    
    % Stage selection
    curr_stage_sfc = 3;
    curr_stage_sp = 3;
    if do_units == 1; curr_stage_sfc = 4; end
    [wrkspc_buffer, group_d, groupm_d ] = Fig_4c2_Bndry_redo(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,opts_PM3Dcs,opts_PM3Dsp,plot_on_func)
    
    group = [group_s, group_d];
    groupm = [groupm_s, groupm_d];

    inds1=[6,5,8,7,2,1,4,3];
    inds2=[3,4,1,2];
    
    if plot_on_spect
        
        opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
        group_pl = group;
        
        % % Sample
        inds = [6,5];
        figure;
        [h1] = plot_matrix3D_custstruct([],group_pl(inds),opts_PM3D,opts_PM3Dcs);
        
        inds = [8,7];
        figure;
        [h1] = plot_matrix3D_custstruct([],group_pl(inds),opts_PM3D,opts_PM3Dcs);
        
        inds = [2,1];
        figure;
        [h1] = plot_matrix3D_custstruct([],group_pl(inds),opts_PM3D,opts_PM3Dcs);
        
        inds = [4,3];
        figure;
        [h1] = plot_matrix3D_custstruct([],group_pl(inds),opts_PM3D,opts_PM3Dcs);
        
        figure; 
        [h1] = plot_matrix3D_custstruct([],groupm(inds2),opts_PM3D,opts_PM3Dcs);
        
        % % Delay
        inds = [6,5]+12;
        figure;
        [h1] = plot_matrix3D_custstruct([],group_pl(inds),opts_PM3D,opts_PM3Dcs);
        
        inds = [8,7]+12;
        figure;
        [h1] = plot_matrix3D_custstruct([],group_pl(inds),opts_PM3D,opts_PM3Dcs);
        
        
        inds = [2,1]+12;
        figure;
        [h1] = plot_matrix3D_custstruct([],group_pl(inds),opts_PM3D,opts_PM3Dcs);
        
        inds = [4,3]+12;
        figure;
        [h1] = plot_matrix3D_custstruct([],group_pl(inds),opts_PM3D,opts_PM3Dcs);
        
        figure; 
        [h1] = plot_matrix3D_custstruct([],groupm(inds2+6),opts_PM3D,opts_PM3Dcs);

    end
    
    
    % Test bargraph
    if plot_on_bargraph
        % All data
            % 50% sample vs delay
        fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group([inds1(1:4) inds1(1:4)+12]),opts_PSC);
            % 60% sample vs delay
        fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group([inds1(5:8) inds1(5:8)+12]),opts_PSC);
        
        % % % % % Rel vs irrel
            % 50% sample and delay
        fign; [h1, h, p, out.PSC] = plot_stats_custstruct(groupm([inds2(1:2) inds2(1:2)+6]),opts_PSC);
            % 60% sample and delay
        fign; [h1, h, p, out.PSC] = plot_stats_custstruct(groupm([inds2(3:4) inds2(3:4)+6]),opts_PSC);
        
        % % % % % Sample vs delay
            % 60% sample vs delay
        fign; [h1, h, p, out.PSC] = plot_stats_custstruct(groupm([inds2(1) inds2(1)+6]),opts_PSC);
            % 60% sample vs delay
        fign; [h1, h, p, out.PSC] = plot_stats_custstruct(groupm([inds2(3) inds2(3)+6]),opts_PSC);
    end    


%% Fig 4D - Boundary Ctgsetli mode 5
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig4d';

plot_on_spect = 1;
plot_on_bargraph = 1;
plot_on_func = 0;

show_only_diff = 1;             % 1-shows only percent changes
                                % 0-shows original traces
                                
% Stage selection
curr_stage_sfc = 2;
curr_stage_sp = 2;

do_units = 0;

% SFC mode
if do_units==0
    sfc_mode = 22.4514111;
    perm_mode = 22.4514111;
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
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    
    [wrkspc_buffer, group, groupm ] = Fig_4d_Bndry5(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,opts_PM3Dcs,opts_PM3Dsp,plot_on_func,do_units)

    if plot_on_spect
        opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};

        if ~show_only_diff
            inds = [2,8,1];
            figure;
            [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);

            inds = [4,6,3];
            figure;
            [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
        else
            inds = [1,2,4,3];
            figure;
            [h1] = plot_matrix3D_custstruct([],groupm(inds),opts_PM3D,opts_PM3Dcs);
        end

    end
    
    ind1=[2,1,8,7,4,3,6,5];
    ind2=[1,2,4,3];

    % Test bargraph
    if plot_on_bargraph
        if ~show_only_diff
            fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group(ind1),opts_PSC);
        end
%         hmask=false(4);
%         hmask(1,2)=true;hmask(3,4)=true;hmask(1,3)=true;hmask(2,4)=true;
%         opts_PSC.hmask = hmask;
        fign; [h1, h, p, out.PSC] = plot_stats_custstruct(groupm(ind2),opts_PSC);
    end
%% Fig 4E - Boundary Ctgsetli mode 5
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig4e';

plot_on_spect = 1;
plot_on_bargraph = 1;
plot_on_func = 0;

show_only_diff = 0;             % 1-shows only percent changes
                                % 0-shows original traces

do_units = 0;

% SFC mode
if do_units==0
    sfc_mode = 22.451411103;
    perm_mode = 22.451411103;
%     sfc_mode = 22.451511113;    % Partial coherence
%     perm_mode = 22.451511113;
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
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    
    % Stage selection
    curr_stage_sfc = 2;
    curr_stage_sp = 2;
    [wrkspc_buffer, group_s, groupm_s ] = Fig_4d_Bndry5(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,opts_PM3Dcs,opts_PM3Dsp,plot_on_func,do_units)
    
    % Stage selection
    curr_stage_sfc = 3;
    curr_stage_sp = 3;
    [wrkspc_buffer, group_d, groupm_d ] = Fig_4d_Bndry5(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,opts_PM3Dcs,opts_PM3Dsp,plot_on_func,do_units)
    
    group = [group_s, group_d];
    groupm = [groupm_s, groupm_d];
    
    
    if plot_on_spect
        opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};

        if ~show_only_diff
            % % Sample
            inds = [2,1];
            figure;
            [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);

            inds = [4,3];
            figure;
            [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
            
            % % Delay
            inds = [2,1]+8;
            figure;
            [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);

            inds = [4,3]+8;
            figure;
            [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
        else
            % % Sample
            inds = [1,2];
            figure;
            [h1] = plot_matrix3D_custstruct([],groupm(inds),opts_PM3D,opts_PM3Dcs);
            
            % % Delay
            inds = [1,2]+4;
            figure;
            [h1] = plot_matrix3D_custstruct([],groupm(inds),opts_PM3D,opts_PM3Dcs);
        end
    end
    
    ind1=[2,1,2,1,4,3,4,3]; ind1([3,4,7,8]) = ind1([3,4,7,8])+8;
    ind2=[1,2,1,2]; ind2([3,4]) = ind2([3,4]) + 4;

    % Test bargraph
    if plot_on_bargraph
        if ~show_only_diff
            fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group(ind1),opts_PSC);
        else
        hmask=false(4);
        hmask(1,2)=true;hmask(3,4)=true;hmask(1,3)=true;hmask(2,4)=true;
        opts_PSC.hmask = hmask;
            fign; [h1, h, p, out.PSC] = plot_stats_custstruct(groupm(ind2),opts_PSC);
        end
    end


%% Fig 4F - Boundary Ctgsetli mode 5; pFFC vs FFC
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig4f';

plot_on_spect = 1;
plot_on_bargraph = 1;
plot_on_func = 0;

show_only_diff = 1;             % 1-shows only percent changes
                                % 0-shows original traces

do_units = 0;

% SFC mode
if do_units==0
    sfc_mode1 = 22.451411103;   % FFC
    perm_mode1 = 22.451411103;
    sfc_mode2 = 22.451511113;    % Partial coherence
    perm_mode2 = 22.451511113;
% sfc_mode = 2.3015111;
% perm_mode = 2.3015111;
% perm_mode = 5;
% sfc_mode = 41.6514111;
% perm_mode = 41.6514111;
elseif do_units == 1
    sfc_mode1 = 22.451411103;   % FFC
    perm_mode1 = 22.451411103;
    sfc_mode2 = 52.7502010;     % Units
    perm_mode2 = 52.7500010;
elseif do_units == 2
    sfc_mode1 = 22.451411103;   % FFC
    perm_mode1 = 22.451411103;
    sfc_mode2 = 52.7500010;      % Units
    perm_mode2 = 52.7500010;
end

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
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    
    % Stage selection
    curr_stage_sfc = 2;
    curr_stage_sp = 2;
    [wrkspc_buffer, group1_s, group1m_s ] = Fig_4d_Bndry5(wrkspc_buffer,sfc_mode1,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode1,curr_stage_sp,freqband_stats_perm,opts_perm,opts_PM3Dcs,opts_PM3Dsp,plot_on_func,do_units)
    [wrkspc_buffer, group2_s, group2m_s ] = Fig_4d_Bndry5(wrkspc_buffer,sfc_mode2,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode2,curr_stage_sp,freqband_stats_perm,opts_perm,opts_PM3Dcs,opts_PM3Dsp,plot_on_func,do_units)
    
    % Stage selection
    curr_stage_sfc = 3;
    curr_stage_sp = 3;
    [wrkspc_buffer, group1_d, group1m_d ] = Fig_4d_Bndry5(wrkspc_buffer,sfc_mode1,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode1,curr_stage_sp,freqband_stats_perm,opts_perm,opts_PM3Dcs,opts_PM3Dsp,plot_on_func,do_units)
    [wrkspc_buffer, group2_d, group2m_d ] = Fig_4d_Bndry5(wrkspc_buffer,sfc_mode2,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode2,curr_stage_sp,freqband_stats_perm,opts_perm,opts_PM3Dcs,opts_PM3Dsp,plot_on_func,do_units)
    
    group1 = [group1_s, group1_d];
    group1m = [group1m_s, group1m_d];
    
    group2 = [group2_s, group2_d];
    group2m = [group2m_s, group2m_d];
    
    
    if plot_on_spect
        opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};

        if ~show_only_diff
            % % Sample
            inds = [2,1];
            figure;
            [h1] = plot_matrix3D_custstruct([],group1(inds),opts_PM3D,opts_PM3Dcs);

            inds = [4,3];
            figure;
            [h1] = plot_matrix3D_custstruct([],group1(inds),opts_PM3D,opts_PM3Dcs);
            
            % % Delay
            inds = [2,1]+8;
            figure;
            [h1] = plot_matrix3D_custstruct([],group1(inds),opts_PM3D,opts_PM3Dcs);

            inds = [4,3]+8;
            figure;
            [h1] = plot_matrix3D_custstruct([],group1(inds),opts_PM3D,opts_PM3Dcs);
            
            % % Sample
            inds = [2,1];
            figure;
            [h1] = plot_matrix3D_custstruct([],group2(inds),opts_PM3D,opts_PM3Dcs);

            inds = [4,3];
            figure;
            [h1] = plot_matrix3D_custstruct([],group2(inds),opts_PM3D,opts_PM3Dcs);
            
            % % Delay
            inds = [2,1]+8;
            figure;
            [h1] = plot_matrix3D_custstruct([],group2(inds),opts_PM3D,opts_PM3Dcs);

            inds = [4,3]+8;
            figure;
            [h1] = plot_matrix3D_custstruct([],group2(inds),opts_PM3D,opts_PM3Dcs);
        else
            % % Sample
            inds = [1,2];
            figure;
            [h1] = plot_matrix3D_custstruct([],group1m(inds),opts_PM3D,opts_PM3Dcs);
            
            % % Delay
            inds = [1,2]+4;
            figure;
            [h1] = plot_matrix3D_custstruct([],group1m(inds),opts_PM3D,opts_PM3Dcs);
            
            % % Sample
            inds = [1,2];
            figure;
            [h1] = plot_matrix3D_custstruct([],group2m(inds),opts_PM3D,opts_PM3Dcs);
            
            % % Delay
            inds = [1,2]+4;
            figure;
            [h1] = plot_matrix3D_custstruct([],group2m(inds),opts_PM3D,opts_PM3Dcs);
        end
    end
    
    
    ind1=[2,1,2,1,4,3,4,3]; ind1([3,4,7,8]) = ind1([3,4,7,8])+8;
    ind2=[1,2,1,2]; ind2([3,4]) = ind2([3,4]) + 4;

    % Test bargraph
    if plot_on_bargraph
        if ~show_only_diff
            fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group1(ind1) group2(ind1)],opts_PSC);
        else
%         hmask=false(4);
%         hmask(1,2)=true;hmask(3,4)=true;hmask(1,3)=true;hmask(2,4)=true;
%         opts_PSC.hmask = hmask;
            fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group1m(1) group2m(1) group1m(2) group2m(2)],opts_PSC);
            fign; [h1, h, p, out.PSC] = plot_stats_custstruct([group1m(1+4) group2m(1+4) group1m(2+4) group2m(2+4)],opts_PSC);
        end
    end

%% Fig 4G - 2xbnd vs 1xbnd (reviewer requested figure)
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig4g';

addpath('figs_batch_5_0_supporting');

merge_rel_irrel = 0;
    show_only_relevant = 0;

% Setup params 
    clear group group0

    data_mode = 22;
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
            if merge_rel_irrel
                s.sfc_mode =  22.461810303;
            else
                s.sfc_mode =  22.461810302;
            end
            s.perm_mode = s.sfc_mode;
        case 23         % FFC spectrogram
            s.sfc_mode =  23.4014111000;
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
    opts_exclude.excludeL = 0;
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

    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;      % 0-Either; 1-Both; 2-Positive; 3-Negative
    opts_perm.alpha0 = 0.05;
    opts_perm.alpha_bh0 = 0.2;
    opts_perm.alpha_bh0 = 0.05;
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
    s.plot_on_spect = 0;
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
    s.opts_PM3Dsp.show_range_stats = 0;
    s.opts_PM3Dsp.show_range_perm = 0;
    s.opts_PSC.paperfig_mode = paperfig_mode;
    s.opts_PSC.remove_dependent = 0;
    s.opts_PSC.do_stats = 1;
    s.opts_PSC.do_diagtest = 0;
    s.opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    s.opts_PSC.opts_diagtest.threshold = 1.0; 
    s.opts_PDays.do_subplots = 1;
        s.opts_PDays.max_subplots_per_fig = 16;
        

[wrkspc_buffer, out] = Fig4g_bnd2x(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)


%
opts_PSC.remove_dependent = 1;
group = out.group
if merge_rel_irrel
    inds = [1,2,10,3,4,10];
    figure;[h1] = plot_matrix3D_custstruct([],group(inds(1:3)),opts_PM3D,opts_PM3Dcs);
    figure;[h1] = plot_matrix3D_custstruct([],group(inds(4:6)),opts_PM3D,opts_PM3Dcs);

    opts_PSC.hmask = ones(3);
    figure('Position',[421    49   401   757]); [h1, h, p, out.PSC] = plot_stats_custstruct(group(inds(1:3)),opts_PSC);
    figure('Position',[421    49   401   757]); [h1, h, p, out.PSC] = plot_stats_custstruct(group(inds(4:6)),opts_PSC);

else
    % Relevant Sch A & B
    inds = [1,2,10,3,4,12];
    figure;[h1] = plot_matrix3D_custstruct([],group(inds(1:3)),opts_PM3D,opts_PM3Dcs);
    figure;[h1] = plot_matrix3D_custstruct([],group(inds(4:end)),opts_PM3D,opts_PM3Dcs);

    opts_PSC.hmask = ones(3);
    figure('Position',[421    49   401   757]); [h1, h, p, out.PSC] = plot_stats_custstruct(group(inds(1:3)),opts_PSC);
    figure('Position',[421    49   401   757]); [h1, h, p, out.PSC] = plot_stats_custstruct(group(inds(4:end)),opts_PSC);
    
    if ~show_only_relevant
        % Irrelevant Sch A & B
        inds = [7,8,12,5,6,10];
        figure;[h1] = plot_matrix3D_custstruct([],group(inds(1:3)),opts_PM3D,opts_PM3Dcs);
        figure;[h1] = plot_matrix3D_custstruct([],group(inds(4:end)),opts_PM3D,opts_PM3Dcs);

        opts_PSC.hmask = ones(3);
        figure('Position',[421    49   401   757]); [h1, h, p, out.PSC] = plot_stats_custstruct(group(inds(1:3)),opts_PSC);
        figure('Position',[421    49   401   757]); [h1, h, p, out.PSC] = plot_stats_custstruct(group(inds(4:end)),opts_PSC);
    end
end



% opts_PSC.hmask = blkdiag(ones(3),ones(3));
% fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group(inds),opts_PSC);
    
%% Fig 15 - FFC normal vs partial coherence
% Results depending on how run sims. See results plots folder

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);

currfigname = 'Fig15';

% Stage selection
curr_stage_sfc = 2;
curr_stage_sp = 2;

% Set up pls fv (popts)
% Unit exclusion
opts_exclude = Opts_Exclude;
opts_exclude.exclude_clipping = 1;
opts_exclude.exclude_60 = 0;
opts_exclude.exclude_nans = 1;
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 0; 
opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC


% More pls switches
opts_pls = Opts_Pls;
opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
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
opts_pls.target_pls_format = 0; % Convert pls to match this format!
opts_pls.collapse_pls_to_days = 0;


% Permutation test options
opts_perm = Opts_Perm;
opts_perm.do_bh0 = 0;
opts_perm.do_phi = 0;
opts_perm.split_plusminus = 0;

% Plotting options
opts_PSC.paperfig_mode = paperfig_mode;
opts_PSC.remove_dependent = 0;
opts_PSC.do_stats = 1;
opts_PSC.do_diagtest = 0;

% Load Ctg response
% SFC mode

sfc_mode =  22.451511113;
perm_mode = 22.451411103;

% sfc_mode =  22.441511112;
% perm_mode = 22.441411102;

opts_pls.permdat2pls = 1;
opts_pls.perm2pls = 0;
opts_pls.do_diff = 1;
% opts_pls.permdat2pls = 0;
% opts_pls.perm2pls = 1;
% opts_pls.do_diff = 0;

curr_stage_sfc = 2; curr_stage_sp = 2;
[wrkspc_buffer, group, out] = Fig_15_FFC_partial_coherence(wrkspc_buffer, sfc_mode, perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,opts_pls,opts_perm,opts_PSC);

curr_stage_sfc = 3; curr_stage_sp = 3;
[wrkspc_buffer, group, out] = Fig_15_FFC_partial_coherence(wrkspc_buffer, sfc_mode, perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,opts_pls,opts_perm,opts_PSC);

%
sfc_mode =  22.401511111;
perm_mode = 22.401511101;

opts_pls.permdat2pls = 0;
opts_pls.perm2pls = 1;
opts_pls.do_diff = 0;
opts_pls.perm2pls_return_mode = 4;

curr_stage_sfc = 2; curr_stage_sp = 2;
[wrkspc_buffer, group, out] = Fig_15_FFC_partial_coherence(wrkspc_buffer, sfc_mode, perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,opts_pls,opts_perm,opts_PSC);

%
curr_stage_sfc = 3; curr_stage_sp = 3;
[wrkspc_buffer, group, out] = Fig_15_FFC_partial_coherence(wrkspc_buffer, sfc_mode, perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,opts_pls,opts_perm,opts_PSC);



%% Fig 15b - FFC partial coherence

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);

currfigname = 'Fig15b';


% SFC mode
sfc_mode =  22.401511111;
perm_mode = 22.401511111;


% Set up pls fv (popts)
% Unit exclusion
opts_exclude = Opts_Exclude;
opts_exclude.exclude_clipping = 1;
opts_exclude.exclude_60 = 0;
opts_exclude.exclude_nans = 1;
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 0; 
opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC

% More pls switches
opts_pls = Opts_Pls;
opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
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
opts_pls.target_pls_format = 0; % Convert pls to match this format!
opts_pls.collapse_pls_to_days = 0;


% Permutation test options
opts_perm = Opts_Perm;
opts_perm.do_bh0 = 1;
opts_perm.do_phi = 0;
opts_perm.split_plusminus = 0;
opts_perm.alpha0 = 0.05;
opts_perm.alpha_bh0 = 0.2;

% Plotting options
opts_PSC.paperfig_mode = 1;
opts_PSC.remove_dependent = 0;
opts_PSC.do_stats = 1;
opts_PSC.do_diagtest = 0;


% Stage selection
curr_stage_sfc = 2;
curr_stage_sp = 2;

opts_pls.perm2pls_return_mode = 4;
opts_PSC.do_stats = 1;
opts_PSC.do_diagtest = 0;
[wrkspc_buffer, group, out] = Fig_15b_FFC_partial_coherence_all(wrkspc_buffer, sfc_mode, perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,opts_pls,opts_perm,opts_PSC);


opts_pls.perm2pls_return_mode = 1;
opts_PSC.do_stats = 0;
opts_PSC.do_diagtest = 1;
[wrkspc_buffer, group, out] = Fig_15b_FFC_partial_coherence_all(wrkspc_buffer, sfc_mode, perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,opts_pls,opts_perm,opts_PSC);


% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 3;

opts_pls.perm2pls_return_mode = 4;
opts_PSC.do_stats = 1;
opts_PSC.do_diagtest = 0;
[wrkspc_buffer, group, out] = Fig_15b_FFC_partial_coherence_all(wrkspc_buffer, sfc_mode, perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,opts_pls,opts_perm,opts_PSC);


opts_pls.perm2pls_return_mode = 1;
opts_PSC.do_stats = 0;
opts_PSC.do_diagtest = 1;
[wrkspc_buffer, group, out] = Fig_15b_FFC_partial_coherence_all(wrkspc_buffer, sfc_mode, perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,opts_pls,opts_perm,opts_PSC);


%% Figure 16a - More detailed code for playing around with of 60% and 100% responses
% Tweaked version of Fig13a

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig16a';


plot_on_spect = 1;
plot_on_bargraph = 0;
plot_on_func = 0;

bri_mode = 0;   % 0-Just plot traces; 1-Difference; 2-Contrast index

do_units=5;             % 0-FFC, units sp;
                        % 1-Units time series, FFC bndry SP
                        % 2-FR traces ctgmode=7, ctgmode=0 SP
                        % 3-FR traces ctgmode=7, ctgmode=7 SP
                        % 4-FR traces ctgmode=7, ctgmode=7 SP - only show "preferred" traces - demonstrate 60%'s are "late bloomers"
                        % 5-all units




switch do_units
    case 0

            groupmode = 0; 
                    % 0 - all cells
                    % 1 - grouped - works for ctgsetli mode 0, 5
                    % 2 - grouped - 60% vs 100%'s for ctgsetli mode 7
                    % 3 - grouped - 60% vs 100%'s for ctgsetli mode 7


            % SFC mode
            sfc_mode =  22.451511103;
%             perm_mode = 22.441511100;
            % sfc_mode = 41.6515111;
            % % perm_mode = 41.6013111;
            % sfc_mode = 52.770001001;          % For looking at unit response
%             sfc_mode = 52.770201001;          % For looking at unit response
            perm_mode = 52.770001001;

            % Stage selection
            curr_stage_sp = 2;
            curr_stage_sfc = 2;
        

    case 1

        groupmode = 1; 

        % SFC mode
%         sfc_mode =  22.451511103;
        perm_mode = 22.451511100;
%         perm_mode = 22.451511103;
        % sfc_mode = 41.6515111;
        % % perm_mode = 41.6013111;
        % sfc_mode = 52.770001001;          % For looking at unit response
        sfc_mode = 52.770201001;          % For looking at unit response
%         perm_mode = 52.770001001;

        % Stage selection
        curr_stage_sp = 2;
        curr_stage_sfc = 4;
            
            
    case 2
        
        groupmode = 1; 
        % SFC mode
        sfc_mode = 52.770201001;          % For looking at unit response
        perm_mode = 52.700001001;

        % Stage selection
        curr_stage_sp = 2;
        curr_stage_sfc = 4;
        
    case 3
        
        groupmode = 3; 
        % SFC mode
        sfc_mode = 52.770201001;          % For looking at unit response
        perm_mode = 52.770001001;

        % Stage selection
        curr_stage_sp = 2;
        curr_stage_sfc = 4;
        
    case 4
        
        groupmode = 4; 
        % SFC mode
        sfc_mode = 52.770201001;          % For looking at unit response
        perm_mode = 52.770001001;

        % Stage selection
        curr_stage_sp = 2;
        curr_stage_sfc = 4;
            
    case 5
        
        groupmode = 5; 
        % SFC mode
        sfc_mode = 52.770201001;          % For looking at unit response
        perm_mode = 52.770001001;

        % Stage selection
        curr_stage_sp = 2;
        curr_stage_sfc = 4;
        
        bri_mode = 0;
            
        
end


    % Pls switches
    freqband_stats = [16 20];
    freqband_stats_perm = [16 20];


    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 0;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 1;
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;  % Only take negatives
    opts_perm.do_quantiles_mode = 0;
            opts_perm.chosen_quantile = 0.1;
            opts_perm.upper_quantile = 0;
            
    if do_units == 1
        % Permutation test options
        opts_perm = Opts_Perm;
        opts_perm.do_bh0 = 1;
        opts_perm.do_phi = 0;
        opts_perm.split_plusminus = 3;  % Only take negatives
        opts_perm.do_quantiles_mode = 0;
                opts_perm.chosen_quantile = 0.1;
                opts_perm.upper_quantile = 0;
    end
    
    
% Load pls
[wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
pls = out_pls.pls;
pls_stats = out_pls.pls_stats;
abscissa = out_pls.abscissa;
bad_any = out_pls.bad_any;
funames1D = out_pls.funames1D;
mypairs = out_pls.mypairs;
group0 = out_pls.group;
clear out_pls


% Load sp's
[wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm,opts_exclude);
sp = out_perm.sig_cells;

% Load bads perm
[bad_any_perm] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);

% Map sp's as needed
[sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any);


% Create group template
mycrit = [2*ones(1,size(sp,2))];
grt = group0(1);
grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
grt.xdata = group0(1).xdata;

switch groupmode 
    
    case 1

        % Run a simple test
        clear group
        i=0;
        i=i+1; group(i)=grt; group(i).criteria([1])=[1]; group(i).ctgs=5;   % Ctg1-2 boundary nboundary ctgs
        i=i+1; group(i)=grt; group(i).criteria([1])=[1]; group(i).ctgs=1;   % Ctg1-2 boundary noundary ctgs
        i=i+1; group(i)=grt; group(i).criteria(1)=[0]; group(i).ctgs=5;   % Ctg1-2 non-deciders
        i=i+1; group(i)=grt; group(i).criteria(1)=[0]; group(i).ctgs=1;   % Ctg1-2 non-deciders

        i=i+1; group(i)=grt; group(i).criteria([2])=[1]; group(i).ctgs=6;   % Ctg1-2 deciders
        i=i+1; group(i)=grt; group(i).criteria([2])=[1]; group(i).ctgs=2;   % Ctg1-2 deciders
        i=i+1; group(i)=grt; group(i).criteria(2)=[0]; group(i).ctgs=6;   % Ctg1-2 non-deciders
        i=i+1; group(i)=grt; group(i).criteria(2)=[0]; group(i).ctgs=2;   % Ctg1-2 non-deciders

    
    case 0
        
        % Run a simple test
        clear group
        i=0;
        i=i+1; group(i)=grt; group(i).criteria([1,5])=[1,2]; group(i).ctgs=1;   % Ctg1-2 boundary nboundary ctgs
        i=i+1; group(i)=grt; group(i).criteria([1,5])=[2,1]; group(i).ctgs=1;   % Ctg1-2 boundary noundary ctgs

        i=i+1; group(i)=grt; group(i).criteria([2,6])=[1,2]; group(i).ctgs=2;   % Ctg1-2 deciders
        i=i+1; group(i)=grt; group(i).criteria([2,6])=[2,1]; group(i).ctgs=2;   % Ctg1-2 deciders

        
    case 3

        % Run a simple test
        clear group
        i=0;
        i=i+1; group(i)=grt; group(i).criteria(5)=[1]; group(i).ctgs=5;   % Ctg1-2 boundary nboundary ctgs
        i=i+1; group(i)=grt; group(i).criteria(5)=[1]; group(i).ctgs=1;   % Ctg1-2 boundary noundary ctgs
        i=i+1; group(i)=grt; group(i).criteria(1)=[1]; group(i).ctgs=5;   % Ctg1-2 non-deciders
        i=i+1; group(i)=grt; group(i).criteria(1)=[1]; group(i).ctgs=1;   % Ctg1-2 non-deciders

        i=i+1; group(i)=grt; group(i).criteria(6)=[1]; group(i).ctgs=6;   % Ctg1-2 deciders
        i=i+1; group(i)=grt; group(i).criteria(6)=[1]; group(i).ctgs=2;   % Ctg1-2 deciders
        i=i+1; group(i)=grt; group(i).criteria(2)=[1]; group(i).ctgs=6;   % Ctg1-2 non-deciders
        i=i+1; group(i)=grt; group(i).criteria(2)=[1]; group(i).ctgs=2;   % Ctg1-2 non-deciders
        
    case 4
        

        % Run a simple test
        clear group
        i=0;
        i=i+1; group(i)=grt; group(i).criteria([1,5])=[1,2]; group(i).ctgs=5;   % Ctg1-2 boundary nboundary ctgs
        i=i+1; group(i)=grt; group(i).criteria([1,5])=[1,2]; group(i).ctgs=1;   % Ctg1-2 boundary noundary ctgs
        i=i+1; group(i)=grt; group(i).criteria([1,5])=[0,1]; group(i).ctgs=5;   % Ctg1-2 boundary nboundary ctgs
        i=i+1; group(i)=grt; group(i).criteria([1,5])=[0,1]; group(i).ctgs=5;   % Ctg1-2 boundary noundary ctgs

        i=i+1; group(i)=grt; group(i).criteria([2,6])=[1,2]; group(i).ctgs=6;   % Ctg1-2 deciders
        i=i+1; group(i)=grt; group(i).criteria([2,6])=[1,2]; group(i).ctgs=2;   % Ctg1-2 deciders
        i=i+1; group(i)=grt; group(i).criteria([2,6])=[0,1]; group(i).ctgs=6;   % Ctg1-2 deciders
        i=i+1; group(i)=grt; group(i).criteria([2,6])=[0,1]; group(i).ctgs=2;   % Ctg1-2 deciders



        % Run a simple test
        clear group
        i=0;
        i=i+1; group(i)=grt; group(i).criteria([1,5])=[2,1]; group(i).ctgs=5;   % Ctg1-2 boundary nboundary ctgs
        i=i+1; group(i)=grt; group(i).criteria([1,5])=[1,2]; group(i).ctgs=1;   % Ctg1-2 boundary noundary ctgs

        i=i+1; group(i)=grt; group(i).criteria([2,6])=[2,1]; group(i).ctgs=6;   % Ctg1-2 deciders
        i=i+1; group(i)=grt; group(i).criteria([2,6])=[1,2]; group(i).ctgs=2;   % Ctg1-2 deciders
        
        
    case 5


        % Run a simple test
        clear group
        i=0;
        i=i+1; group(i)=grt; group(i).criteria(1)=[2]; group(i).ctgs=5;   % Ctg1-2 boundary nboundary ctgs
        i=i+1; group(i)=grt; group(i).criteria(1)=[2]; group(i).ctgs=1;   % Ctg1-2 boundary noundary ctgs
        i=i+1; group(i)=grt; group(i).criteria(1)=[2]; group(i).ctgs=5;   % Ctg1-2 non-deciders
        i=i+1; group(i)=grt; group(i).criteria(1)=[2]; group(i).ctgs=1;   % Ctg1-2 non-deciders

        i=i+1; group(i)=grt; group(i).criteria(2)=[2]; group(i).ctgs=6;   % Ctg1-2 deciders
        i=i+1; group(i)=grt; group(i).criteria(2)=[2]; group(i).ctgs=2;   % Ctg1-2 deciders
        i=i+1; group(i)=grt; group(i).criteria(2)=[2]; group(i).ctgs=6;   % Ctg1-2 non-deciders
        i=i+1; group(i)=grt; group(i).criteria(2)=[2]; group(i).ctgs=2;   % Ctg1-2 non-deciders
        
        
        
        

end

% Calculate legend entries
group = group.query_legend(group0);

if bri_mode == 2
    % Calculate boundary resolvability index
    group = grouppairs_merge(group);
    pls_bnd = pls(:,:,1:4); pls_nbnd = pls(:,:,5:8);
    pls_bri = (pls_bnd - pls_nbnd) ./ (pls_nbnd + pls_bnd);
    pls = pls_bri;
    i=0;
    i=i+1; group(i).ctgs = 1; group(i).legend = 'Bnd BRI Cat vs Dog';
    i=i+1; group(i).ctgs = 1; group(i).legend = 'Non-Bnd BRI Cat vs Dog';
    i=i+1; group(i).ctgs = 2; group(i).legend = 'Bnd BRI Goc vs Tad';
    i=i+1; group(i).ctgs = 2; group(i).legend = 'Non-Bnd BRI Goc vs Tad';
end

% Load data into groups
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D);


% TAke difference
if bri_mode == 1
    group = grouppairs_merge(group);
end

% Test Plot groups, at last

if plot_on_spect
    figure;
    inds = 1:round(length(group)/2);
    [h1] = plot_matrix3D_custstruct([],group(inds));
    
    figure;
    inds = round(length(group)/2)+1:length(group);
    [h1] = plot_matrix3D_custstruct([],group(inds));
    
end

% 
opts_PSC.remove_dependent = 0;
if plot_on_bargraph
    
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
end



%% Figure 16b - As 16a, but groups criteria reflects both two SP's - FFC and unit activity

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig16b';


plot_on_spect = 1;
plot_on_bargraph = 0;
plot_on_func = 0;

bri_mode = 0;   % 0-Just plot traces; 1-Difference; 2-Contrast index

do_units = 1            % 0-FFC, units sp;
                        % 1-Units time series, FFC bndry SP
                        % 2-FR traces ctgmode=7, ctgmode=0 SP
                        % 3-FR traces ctgmode=7, ctgmode=7 SP
                        % 4-FR traces ctgmode=7, ctgmode=7 SP - only show "preferred" traces - demonstrate 60%'s are "late bloomers"
                        % 5-all units

switch do_units
    case 1

        groupmode = 1; 

        % SFC mode
%         sfc_mode =  22.451511103;
        perm_mode2 = 22.451511100;
%         perm_mode = 22.451511103;
        % sfc_mode = 41.6515111;
        % % perm_mode = 41.6013111;
        % sfc_mode = 52.770001001;          % For looking at unit response
        sfc_mode = 52.770201001;          % For looking at unit response
        perm_mode = 52.770001001;

        % Stage selection
        curr_stage_sp = 2;
        curr_stage_sp2 = 2;
        curr_stage_sfc = 4;

end

    % Pls switches
    freqband_stats = [16 20];
    freqband_stats_perm = [16 20];


    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 0;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 1;
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    opts_perm.do_quantiles_mode = 0;
            opts_perm.chosen_quantile = 0.1;
            opts_perm.upper_quantile = 0;
            
    warning('do single animals');
    opts_perm2 = opts_perm;
    opts_perm2.split_plusminus = 3;   % Only take negatives
    opts_perm2.do_quantiles_mode = 1;
    
% Load pls
[wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
pls = out_pls.pls;
pls_stats = out_pls.pls_stats;
abscissa = out_pls.abscissa;
bad_any = out_pls.bad_any;
funames1D = out_pls.funames1D;
mypairs = out_pls.mypairs;
group0 = out_pls.group;
clear out_pls


% Load sp's
[wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm,opts_exclude);
sp1 = out_perm.sig_cells;

% Load bads perm
[bad_any_perm] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);

% Map sp's as needed
[sp1] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp1,wrkspc_buffer.currmd.md,bad_any_perm,bad_any);

% Load sp2's
[wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode2,curr_stage_sp2,freqband_stats_perm,bad_any,opts_perm2,opts_exclude);
sp2 = out_perm.sig_cells;

% Load bads perm
[bad_any_perm2] = load_bads(perm_mode2,curr_stage_sp2,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);

% Map sp's as needed
[sp2] = map_sp(perm_mode2, sfc_mode,out_perm.mypairs,mypairs,sp2,wrkspc_buffer.currmd.md,bad_any_perm,bad_any);


% Merge SP's
sp = [sp1 sp2];



% Create group template
mycrit =  2*ones(1,size(sp1,2));
mycrit2 = 2*ones(1,size(sp2,2));

grt = group0(1);
grt.criteria = mycrit; grt.criteria_alt = mycrit2; grt.criteria_sfc = []; grt.ctgs = 1;
grt.xdata = group0(1).xdata;

switch groupmode 
    
    case 1

        % Run a simple test
        clear group
        i=0;
        % Cat/Dog
        i=i+1; group(i)=grt; group(i).criteria([5])=[1]; group(i).criteria_alt([1])=[1]; group(i).ctgs=5;
        i=i+1; group(i)=grt; group(i).criteria([5])=[1]; group(i).criteria_alt([1])=[0]; group(i).ctgs=5;
        i=i+1; group(i)=grt; group(i).criteria([1])=[1]; group(i).criteria_alt([1])=[1]; group(i).ctgs=1;
        i=i+1; group(i)=grt; group(i).criteria([1])=[1]; group(i).criteria_alt([1])=[0]; group(i).ctgs=1;
        
        % Goc/Tad
        i=i+1; group(i)=grt; group(i).criteria([6])=[1]; group(i).criteria_alt([2])=[1]; group(i).ctgs=6;
        i=i+1; group(i)=grt; group(i).criteria([6])=[1]; group(i).criteria_alt([2])=[0]; group(i).ctgs=6;
        i=i+1; group(i)=grt; group(i).criteria([2])=[1]; group(i).criteria_alt([2])=[1]; group(i).ctgs=2;
        i=i+1; group(i)=grt; group(i).criteria([2])=[1]; group(i).criteria_alt([2])=[0]; group(i).ctgs=2;
        
end

% Calculate legend entries
group = group.query_legend(group0);

if bri_mode == 2
    % Calculate boundary resolvability index
    group = grouppairs_merge(group);
    pls_bnd = pls(:,:,1:4); pls_nbnd = pls(:,:,5:8);
    pls_bri = (pls_bnd - pls_nbnd) ./ (pls_nbnd + pls_bnd);
    pls = pls_bri;
    i=0;
    i=i+1; group(i).ctgs = 1; group(i).legend = 'Bnd BRI Cat vs Dog';
    i=i+1; group(i).ctgs = 1; group(i).legend = 'Non-Bnd BRI Cat vs Dog';
    i=i+1; group(i).ctgs = 2; group(i).legend = 'Bnd BRI Goc vs Tad';
    i=i+1; group(i).ctgs = 2; group(i).legend = 'Non-Bnd BRI Goc vs Tad';
end

% Load data into groups
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D);


% TAke difference
if bri_mode == 1
    group = grouppairs_merge(group);
end

% Test Plot groups, at last

if plot_on_spect
    figure;
    inds = 1:round(length(group)/2);
    [h1] = plot_matrix3D_custstruct([],group(inds));
    
    figure;
    inds = round(length(group)/2)+1:length(group);
    [h1] = plot_matrix3D_custstruct([],group(inds));
    
end

% 
opts_PSC.remove_dependent = 0;
if plot_on_bargraph
    
    i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
end

%% Fig 17 - Category vs Bndry Information Index for FFC Normal vs Partial

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);

plot_on_spect=1;
plot_on_bargraph=1;

currfigname = 'Fig17';

% SFC mode FFC
sfc_mode1 =  22.401411101;
sfc_mode2 =  22.451411103;
% perm_mode = 22.401511111;


% SFC mode pFFC
psfc_mode1 =  22.401511111;  % pFFC
psfc_mode2 =  22.451511113;
% perm_mode = 22.401511111;
psfc_mode1 =  52.7000010;    % Uncomment for units
psfc_mode2 =  52.7500010;


% Stage selection
curr_stage_sfc1 = 2;
curr_stage_sfc2 = 2;

    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 1;
    opts_exclude.excludeL = 0;
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC

    warning('To do - implement Buschmann figure.');
    warning('Add horizontal bar in perm figures to indicate value expected by chance');

    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
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
    opts_pls.target_pls_format = 0; % Convert pls to match this format!
    opts_pls.collapse_pls_to_days = 0;
        

    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    opts_perm.alpha0 = 0.05;
    opts_perm.alpha_bh0 = 0.2;
    
    
    [wrkspc_buffer, group_FFC] = Fig_17_FFC_partial_coherence_index(wrkspc_buffer, sfc_mode1, sfc_mode2,curr_stage_sfc1, curr_stage_sfc2,opts_exclude,opts_pls)
    
    [wrkspc_buffer, group_pFFC] = Fig_17_FFC_partial_coherence_index(wrkspc_buffer, psfc_mode1, psfc_mode2,curr_stage_sfc1, curr_stage_sfc2,opts_exclude,opts_pls)
    
    group = [group_FFC(1), group_pFFC(1), group_FFC(2), group_pFFC(2)]; 
    

   
% Test Plot groups, at last
opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
if plot_on_spect
    
    N=length(group);
        inds = 1:N;

    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);

end


% Test bargraph
if plot_on_bargraph
    fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
end



%% Figure 17b - FFC Information Index
% Uses ctgsetli mode 9; based on Fig9b

% Set up default parameters
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig17b';


plot_on_spect=1;
plot_on_bargraph=1;
plot_on_func = 0;

% Setup unitsCI_vs_electBNDRY if we're using it
do_unitsCI_vs_electBNDRY = 0;


% Pls switches
freqband_stats = [16 20];
freqband_stats_perm = [16 20];
% freqband_stats = [10 12];
% freqband_stats_perm = [10 12];

% Set up pls fv (popts)
% Unit exclusion
opts_exclude = Opts_Exclude;
opts_exclude.exclude_clipping = 1;
opts_exclude.exclude_60 = 0;
opts_exclude.exclude_nans = 1;
opts_exclude.excludeL = 0; 
opts_exclude.excludeO = 0; 
opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC


% More pls switches
opts_pls = Opts_Pls;
opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
    opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
    opts_pls.perm2pls_dophi = 0;
    opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
    opts_pls.swap_pls = [];
    %opts_pls.swap_pls = [2,4;5,7];      
    opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
opts_pls.collapse_pls_to_days = 0;

% Permutation test options
opts_perm = Opts_Perm;
opts_perm.do_bh0 = 1;
opts_perm.do_phi = 0;
opts_perm.split_plusminus = 0;


% Stage selection
curr_stage_sp = 2;
curr_stage_sfc = 2;

% SFC mode1 - FFC
sfc_mode = 22.49141030;
perm_mode_sort = 22.401411100;
perm_mode_sp = 22.401411100;
[wrkspc_buffer, gall1] = Fig_17b_Ctgsetli9_Bndry_Information_Index(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode_sort,perm_mode_sp,curr_stage_sp,freqband_stats_perm,opts_perm,do_unitsCI_vs_electBNDRY,plot_on_func)
%figure; [h1] = plot_matrix3D_custstruct([],gall1,opts_PM3D,opts_PM3Dcs);


% SFC mode2 - partialFFC
sfc_mode = 22.491510310;                 % Partial coherence all cells
perm_mode_sort = 22.401511110;
perm_mode_sp =   22.401511110;
% sfc_mode = 52.790000001;                 % Uncomment for units
% perm_mode_sort = 52.700001001;
% perm_mode_sp =   52.700001001;

[wrkspc_buffer, gall2] = Fig_17b_Ctgsetli9_Bndry_Information_Index(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode_sort,perm_mode_sp,curr_stage_sp,freqband_stats_perm,opts_perm,do_unitsCI_vs_electBNDRY,plot_on_func)
%figure; [h1] = plot_matrix3D_custstruct([],gall2,opts_PM3D,opts_PM3Dcs);

for i = 1:length(gall1); gall1(i).legend = ['FFC ' gall1(i).legend]; end
for i = 1:length(gall2); gall2(i).legend = ['pFFC ' gall2(i).legend]; end

group = [gall1(1) gall2(1) gall1(2) gall2(2)]; 
   
% Test Plot groups, at last
opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
if plot_on_spect
    figure;
    [h1] = plot_matrix3D_custstruct([],group,opts_PM3D,opts_PM3Dcs);
end


% Test bargraph
if plot_on_bargraph
    for i = 1:length(group); group(i).legend = ' '; end
    opts_PSC.do_stats = 1;
    fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
end

%% Fig 18a - Electrode distances
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig18a';

addpath('figs_batch_5_0_supporting');


% Setup params 
    clear group group0

    data_mode = 22;
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
            s.sfc_mode =  22.401411101;
            s.perm_mode = s.sfc_mode;
%             s.sfc_mode =  22.401511111;     % Partial FFC
%             s.perm_mode = s.sfc_mode;
%             s.sfc_mode =  22.451411103;     % FFC boundary vs non-boundary
%             s.perm_mode = s.sfc_mode;
        case 23         % FFC spectrogram
            s.sfc_mode =  23.481810103;
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
            s.sfc_mode  = 52.700201001;
            s.perm_mode = 52.700001001;
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
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % Pls switches
    s.freqband_stats = [16 20];
    s.freqband_perm = [16 20];
%     s.freqband_stats = [10 12];
%     s.freqband_perm = [10 12];
    s.timeband_stats = [.6]; s.timeband_perm = [.6];
    s.tf_label_stats = 'Default'; s.tf_label_perm = 'Default';
    [tf_avail] = get_freqband_timeband(s.sfc_mode,opts_exclude); s.tf_avail = tf_avail;
    
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

    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;      % 0-Either; 1-Both; 2-Positive; 3-Negative
    opts_perm.alpha0 = 0.05;
    opts_perm.alpha_bh0 = 0.2;
%     opts_perm.alpha_bh0 = 0.05;
    opts_perm.do_quantiles_mode = 0;
        opts_perm.chosen_quantile = .15;
        opts_perm.upper_quantile = 0;
    
    % Map sp parameters
    s.sp_threshold = 10;
    
    % Group options
    s.do_group_collapse_pls2days = 0;
    
    % Groupmode
    s.groupmode = 1.5;   % 0-Use default grouping (all pairs, enumerate over ctgs);
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
    s.plot_on_electrode_distances = 1;
    
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
    s.opts_PM3Dsp.show_range_perm = 0;
    s.opts_PSC.paperfig_mode = paperfig_mode;
    s.opts_PSC.remove_dependent = 0;
    s.opts_PSC.do_stats = 1;
    s.opts_PSC.do_diagtest = 0;
    s.opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    s.opts_PSC.opts_diagtest.threshold = 1.0; 
    s.opts_PDays.do_subplots = 1;
        s.opts_PDays.max_subplots_per_fig = 16;
        
    do_AB_symmetry = 1;
        
% % % % % % % % % % % % % % % % % % Both animals % % % % % % % % % % % % % % 
opts_exclude.excludeL = 0;
opts_exclude.excludeO = 0;
% Stage selection
s.curr_stage_sfc = 2;
s.curr_stage_sp = 2;
[wrkspc_buffer, out] = Fig18a_electrode_distances(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm,do_AB_symmetry);

% Stage selection
s.curr_stage_sfc = 3;
s.curr_stage_sp = 3;
[wrkspc_buffer, out] = Fig18a_electrode_distances(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm,do_AB_symmetry);

% %%
% % % % % % % % % % % % % % % % % % % Animal L% % % % % % % % % % % % % % 
% opts_exclude.excludeL = 0;
% opts_exclude.excludeO = 1;
% % Stage selection
% s.curr_stage_sfc = 2;
% s.curr_stage_sp = 2;
% [wrkspc_buffer, out] = Fig18a_electrode_distances(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm,do_AB_symmetry);
% 
% % Stage selection
% s.curr_stage_sfc = 3;
% s.curr_stage_sp = 3;
% [wrkspc_buffer, out] = Fig18a_electrode_distances(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm,do_AB_symmetry);
% 
% % % % % % % % % % % % % % % % % % % Animal O % % % % % % % % % % % % % % 
% opts_exclude.excludeL = 1;
% opts_exclude.excludeO = 0;
% % Stage selection
% s.curr_stage_sfc = 2;
% s.curr_stage_sp = 2;
% [wrkspc_buffer, out] = Fig18a_electrode_distances(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm,do_AB_symmetry);
% 
% % Stage selection
% s.curr_stage_sfc = 3;
% s.curr_stage_sp = 3;
% [wrkspc_buffer, out] = Fig18a_electrode_distances(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm,do_AB_symmetry);

%% Fig 18b - Electrode locations
clearvars -except wrkspc_buffer fv fv3
vars_pull(fv3);
currfigname = 'Fig18b';

addpath('figs_batch_5_0_supporting');


% Setup params 
    clear group group0

    data_mode = 22;
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
            s.sfc_mode =  22.401411101;
            s.perm_mode = s.sfc_mode;
%             s.sfc_mode =  22.401511111;     % Partial FFC
%             s.perm_mode = s.sfc_mode;
%             s.sfc_mode =  22.451411103;     % FFC boundary vs non-boundary
%             s.perm_mode = s.sfc_mode;
        case 23         % FFC spectrogram
            s.sfc_mode =  23.481810103;
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
            s.sfc_mode  = 52.700201001;
            s.perm_mode = 52.700001001;
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
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % Pls switches
    s.freqband_stats = [16 20];
    s.freqband_perm = [16 20];
%     s.freqband_stats = [10 12];
%     s.freqband_perm = [10 12];
    s.timeband_stats = [.6]; s.timeband_perm = [.6];
    s.tf_label_stats = 'Default'; s.tf_label_perm = 'Default';
    [tf_avail] = get_freqband_timeband(s.sfc_mode,opts_exclude); s.tf_avail = tf_avail;
    
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

    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;      % 0-Either; 1-Both; 2-Positive; 3-Negative
    opts_perm.alpha0 = 0.05;
    opts_perm.alpha_bh0 = 0.2;
%     opts_perm.alpha_bh0 = 0.05;
    opts_perm.do_quantiles_mode = 0;
        opts_perm.chosen_quantile = .15;
        opts_perm.upper_quantile = 0;
    
    % Map sp parameters
    s.sp_threshold = 10;
    
    % Group options
    s.do_group_collapse_pls2days = 0;
    
    % Groupmode
    s.groupmode = 1.5;   % 0-Use default grouping (all pairs, enumerate over ctgs);
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
    s.plot_on_electrode_locations = 1;
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
    s.opts_PM3Dsp.show_range_perm = 0;
    s.opts_PSC.paperfig_mode = paperfig_mode;
    s.opts_PSC.remove_dependent = 0;
    s.opts_PSC.do_stats = 1;
    s.opts_PSC.do_diagtest = 0;
    s.opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    s.opts_PSC.opts_diagtest.threshold = 1.0; 
    s.opts_PDays.do_subplots = 1;
        s.opts_PDays.max_subplots_per_fig = 16;

% Stage selection - Sample
s.curr_stage_sfc = 2;
s.curr_stage_sp = 2;
[wrkspc_buffer, out] = Fig18b_electrode_locations(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)


% Stage selection - Delay
s.curr_stage_sfc = 3;
s.curr_stage_sp = 3;
[wrkspc_buffer, out] = Fig18b_electrode_locations(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)
