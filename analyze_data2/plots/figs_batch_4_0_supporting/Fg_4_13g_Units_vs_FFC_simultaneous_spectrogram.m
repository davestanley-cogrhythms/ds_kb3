function [wrkspc_buffer, group] = Fg_4_13g_Units_vs_FFC_simultaneous_spectrogram(wrkspc_buffer,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,freqband_stats_perm2,opts_exclude,groupmode,do_units)

% Plot switches
plot_on_spect = 1;
plot_on_bargraph = 1;


% Clustering method
use_power = 0;


    if do_units == 0            % FFC
        % SFC mode FFC
        pername = 'Units';
        sfc_mode = 22.401811101;
        perm_mode = 52.700001001;    

    elseif do_units == 1        % Units time series
        % SFC mode units
        pername = 'Ens.';
        sfc_mode = 52.700201001;
% %         perm_mode = 22.401811101;
% %         if use_power; perm_mode = 41.601311101; end
        plot_on_bargraph = 0;
        perm_mode = 23.4018111013;
%         perm_mode = 23.4013111011;

    elseif do_units == 2        % Units by stages
        pername = 'Ens.';
        sfc_mode = 52.700001001;
%         perm_mode = 22.401811101;
%         if use_power; perm_mode = 41.601311101; end
        plot_on_spect = 0;
        perm_mode = 23.4018111013;
%         perm_mode = 23.4013111011;
    end

    if do_units == 1
        curr_stage_sfc=4;
    end
    
    % % Calc timeband and freqband for perm
    if ~opts_exclude.excludeL && opts_exclude.excludeO        % Doing animal L
        % For beta
        [tf_avail] = get_freqband_timeband(perm_mode,opts_exclude); s.tf_avail = tf_avail;
        i=3; s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
        freqband_stats_perm = s.freqband_perm;
        timeband_perm = s.timeband_perm;
        tf_label_perm = s.tf_label_perm;

        % For alpha
        [tf_avail] = get_freqband_timeband(perm_mode,opts_exclude); s.tf_avail = tf_avail;
        i=4; s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
        freqband_stats_perm2 = s.freqband_perm;
        timeband_perm2 = s.timeband_perm;
        tf_label_perm2 = s.tf_label_perm;
        clear s
    elseif opts_exclude.excludeL && ~opts_exclude.excludeO        % Doing animal O
        % For beta
        [tf_avail] = get_freqband_timeband(perm_mode,opts_exclude); s.tf_avail = tf_avail;
        i=10; s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
        freqband_stats_perm = s.freqband_perm;
        timeband_perm = s.timeband_perm;
        tf_label_perm = s.tf_label_perm;

        % For alpha
        [tf_avail] = get_freqband_timeband(perm_mode,opts_exclude); s.tf_avail = tf_avail;
        i=11; s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
        freqband_stats_perm2 = s.freqband_perm;
        timeband_perm2 = s.timeband_perm;
        tf_label_perm2 = s.tf_label_perm;
        clear s
    end

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
    opts_pls.target_pls_format = 0; % Convert pls to match this format!
    opts_pls.collapse_pls_to_days = 0;
        

    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 1;
    opts_perm.alpha0 = 0.05;
    opts_perm.alpha_bh0 = 0.2;
    %opts_perm.alpha_bh0 = 0.05;
    
    % Map sp parameters
    sp_threshold = 20;
    
    % Group options
    do_group_collapse_pls2days = 0;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dcs.stats_mode = 0;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    opts_PSC.do_stats = 1;
    opts_PSC.do_diagtest = 0;
    opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    opts_PSC.opts_diagtest.threshold = 0; 
    
 
%% Plot switches
% plot_on_spect = 1;
plot_on_spectrogram = 0;
plot_on_scatter = 0;
% plot_on_bargraph = 0;
plot_on_histogram = 0;

% groupmode = 1;   % 0-Use default grouping (all pairs, enumerate over ctgs);
                 % 1-Group based on permutation test output for diff.
                 % 2-Group based on permutation test output all.
    swap_mode = 0;

group_do_merge = 0;
    grouppmerge_do_percent = 1;


%% Import data if running as a function; define dependent parameters

running_as_script = ~is_calledby;
if ~running_as_script
    pop_workspace(false);
end

        
%% Load pls
[wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
pls = out_pls.pls;
pls_stats = out_pls.pls_stats;
abscissa = out_pls.abscissa;
abscissa2 = out_pls.abscissa2;
bad_any = out_pls.bad_any;
funames1D = out_pls.funames1D;
mypairs = out_pls.mypairs;
group0 = out_pls.group;
sfc_mode = out_pls.sfc_mode;
clear out_pls

% Invert pls if needed
[~, sfc_subgroups] = decode_sfc_mode(sfc_mode);
[~, ~, ctgsetli_mode] = build_sfcmode(sfc_mode, sfc_subgroups);
if any(ctgsetli_mode == [4,5]) && opts_pls.do_diff
    pls = pls * -1;
    pls_stats = pls_stats * -1;
end


%% Build group

% Group based on perm_mode
if ~exist('group','var')
    
    
    
    % Load sp's
    opts_perm.timeband_perm = timeband_perm;
    opts_perm.tf_label_perm = tf_label_perm;
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm,opts_exclude);
    sp1 = out_perm.sig_cells;
    
    
    opts_perm.timeband_perm = timeband_perm2;
    opts_perm.tf_label_perm = tf_label_perm2;
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm2,bad_any,opts_perm,opts_exclude);
    sp2 = out_perm.sig_cells;
    
    sp = [sp1, sp2];

    % Load bads perm
    [bad_any_perm] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);

    % Map sp's as needed
    [sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any,sp_threshold);

    % Create group template
    mycrit = [2*ones(1,size(sp,2))];
    grt = group0(1);
    grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
    
    switch groupmode
        
        case -1         % Mode for replicating Fg4_13e
            i=0;
            i=i+1; group(i)=grt; group(i).criteria([1:2,11:12])=[0 1 2 2]; group(i).ctgs=1;         % Cat/Dog beta; Cat/Dog alpha
            i=i+1; group(i)=grt; group(i).criteria([1:2,11:12])=[0 0 0 0]; group(i).ctgs=1;
            i=i+1; group(i)=grt; group(i).criteria([1:2,11:12])=[2 2 1 0]; group(i).ctgs=1;
            
            i=0;
            i=i+1;group(i).legend = 'Dog beta';
            i=i+1;group(i).legend = 'No pref';
            i=i+1;group(i).legend = 'Cat alpha';
            
        case -2         % Mode for replicating Fg4_13e
            i=0;
            i=i+1; group(i)=grt; group(i).criteria([3:4,13:14])=[0 1 2 2]; group(i).ctgs=2;
            i=i+1; group(i)=grt; group(i).criteria([3:4,13:14])=[0 0 0 0]; group(i).ctgs=2;
            i=i+1; group(i)=grt; group(i).criteria([3:4,13:14])=[2 2 1 0]; group(i).ctgs=2;

            i=0;
            i=i+1;group(i).legend = 'Thin beta';
            i=i+1;group(i).legend = 'No pref';
            i=i+1;group(i).legend = 'Fat alpha';
        
        
        case 1
            i=0;
            i=i+1; group(i)=grt; group(i).criteria([1:2,11:12])=[0 1 0 0]; group(i).ctgs=1;         % Cat/Dog beta; Cat/Dog alpha
            i=i+1; group(i)=grt; group(i).criteria([1:2,11:12])=[0 0 0 0]; group(i).ctgs=1;
            i=i+1; group(i)=grt; group(i).criteria([1:2,11:12])=[0 0 1 0]; group(i).ctgs=1;
            
            %group = group.query_legend(group0);
            i=0;
            i=i+1;group(i).legend = 'Dog beta';
            i=i+1;group(i).legend = 'No pref';
            i=i+1;group(i).legend = 'Cat alpha';
            
        case 2
            i=0;
            i=i+1; group(i)=grt; group(i).criteria([3:4,13:14])=[0 1 0 0]; group(i).ctgs=2;
            i=i+1; group(i)=grt; group(i).criteria([3:4,13:14])=[0 0 0 0]; group(i).ctgs=2;
            i=i+1; group(i)=grt; group(i).criteria([3:4,13:14])=[0 0 1 0]; group(i).ctgs=2;

            %group = group.query_legend(group0);
            i=0;
            i=i+1;group(i).legend = 'Thin beta';
            i=i+1;group(i).legend = 'No pref';
            i=i+1;group(i).legend = 'Fat alpha';
% % % % % % % % % % %  These (cases 3-6) are for comparing the SAME UNITS under
% two different conditions (1) biasd response under preferred scheme
%                          (2) biased response under non-preferred scheme
% % % % % % I tested this for a variety of parameters and nothing worked. % % % 
        case 3
            i=0;
            i=i+1; group(i)=grt; group(i).criteria([1:2,11:12])=[0 1 2 2]; group(i).ctgs=1;         % Cat/Dog beta; Cat/Dog alpha
            i=i+1; group(i)=grt; group(i).criteria([1:2,11:12])=[0 1 2 2]; group(i).ctgs=2;         % Cat/Dog beta; Cat/Dog alpha
            
            i=0;
            i=i+1;group(i).legend = 'Dog beta:Cat/Dog';
            i=i+1;group(i).legend = 'Dog beta:Fat/Thin';
            
        case 4
            i=0;
            i=i+1; group(i)=grt; group(i).criteria([1:2,11:12])=[2 2 1 0]; group(i).ctgs=1;         % Cat/Dog beta; Cat/Dog alpha
            i=i+1; group(i)=grt; group(i).criteria([1:2,11:12])=[2 2 1 0]; group(i).ctgs=2;         % Cat/Dog beta; Cat/Dog alpha
            
            i=0;
            i=i+1;group(i).legend = 'Cat alpha:Cat/Dog';
            i=i+1;group(i).legend = 'Cat alpha:Fat/Thin';
            
        case 5
            i=0;
            i=i+1; group(i)=grt; group(i).criteria([3:4,13:14])=[0 1 2 2]; group(i).ctgs=2;         % Cat/Dog beta; Cat/Dog alpha
            i=i+1; group(i)=grt; group(i).criteria([3:4,13:14])=[0 1 2 2]; group(i).ctgs=1;         % Cat/Dog beta; Cat/Dog alpha
            
            i=0;
            i=i+1;group(i).legend = 'Thin beta:Fat/thin';
            i=i+1;group(i).legend = 'Thin beta:Cat/Dog';
            
        case 6
            i=0;
            i=i+1; group(i)=grt; group(i).criteria([3:4,13:14])=[2 2 1 0]; group(i).ctgs=2;         % Cat/Dog beta; Cat/Dog alpha
            i=i+1; group(i)=grt; group(i).criteria([3:4,13:14])=[2 2 1 0]; group(i).ctgs=1;         % Cat/Dog beta; Cat/Dog alpha
            
            i=0;
            i=i+1;group(i).legend = 'Fat alpha:Fat/Thin';
            i=i+1;group(i).legend = 'Fat alpha:Cat/Dog';
    end


end

%% Load data into groups
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);

if do_group_collapse_pls2days
    group = group_collapse_pls2days(group);
end

if group_do_merge
    % Flip if ctgsetli mode = 4 or 5.
    temp=1:length(group);
    N=length(group);
    temp=flipud(reshape(1:N,2,floor(N/2)));
    temp=temp(:);
    group = grouppairs_merge(group(temp),opts_pls.perm2pls_dophi,grouppmerge_do_percent);
end

%% Test Plot groups, at last
opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
if plot_on_spect
    
    do_split = 0;

    N=length(group);
    
    inds = 1:N;
    
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    xl = xlim;
    hold on; plot([xl(1), xl(2)],[0,0],'k','LineWidth',2)
    
end

if 0
    %%
    figure;
    inds = [1,2,3,4,9,10];
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D);
%     figure;
%     inds = 3:4;
%     [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D);

end

%% Test spectrogram
if plot_on_spectrogram
    plot_spectrogram_custstruct(group(:),opts_PM3Dsp);
    
end


%% Test bargraph

if plot_on_bargraph
    hmask = true(length(group));
    hmask(logical(eye(size(hmask)))) = false;
    opts_PSC.hmask = hmask;
    opts_PSC.multicomparisons_mode = 2;     % 2-default (holm-bonferroni)
    % opts_PSC.multicomparisons_mode = 0; warning('Multiple comparisons turned off');
    fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
end

%% Plot histogram
if plot_on_histogram
    plot_histplot(group,'plot_type','hist_autogrouped');
end

%% Test Plot scattergroups

if plot_on_scatter

    figure; % subplot(221); 
    indx = 1;
    indy = 5;
    plot_scatterpls_and_groups(pls_stats,group([1,5]),bad_any,indx,indy,group0)
    
    
end

do_pca = 0;
if do_pca
   [coeff,score,latent] = princomp([x(:),y(:)]);
    ax1 = coeff(:,1);% *sqrt(latent(1));
    ax2 = coeff(:,2);% *sqrt(latent(2));
    hold on; plot([0 coeff(1,1)],[0 coeff(2,1)],'r','LineWidth',2);
    hold on; plot([0 coeff(1,2)],[0 coeff(2,2)],'r','LineWidth',2); 
end


% %% Slope test
% if plot_Sch_sens
%     pls_stats2 = pls_stats;
%     pls_stats2(:,group(3).ctgs) = pls_stats(:,group(3).ctgs) ./ pls_stats(:,group(1).ctgs);
%     figure; % subplot(221); 
%     plott_fit(pls_stats2(~bad_any,group(1).ctgs),pls_stats2(~bad_any,group(3).ctgs),'k.');
%     hold on; plot_scattergroups(pls_stats2(~bad_any,group(1).ctgs),pls_stats2(~bad_any,group(3).ctgs),[group(1),group(3)],~bad_any);
%     xlabel('FFC pls(:,1)'); ylabel('FFC pls(:,2)');
% end


%% Package outputs

% Package data if this script is called by a parent function
if is_calledby
    out.wrkspc_buffer = wrkspc_buffer;
    out.group = group;
    out.bad_any = bad_any;
    out.pls = pls;
    out.pls_stats = pls_stats;
    out.abscissa = abscissa;
    out.mypairs = mypairs;
end
end