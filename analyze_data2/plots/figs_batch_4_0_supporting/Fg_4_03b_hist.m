function [wrkspc_buffer, group] = Fg_4_03b_hist(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,opts_exclude,groupmode)


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
        opts_pls.do_diff = 1;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        opts_pls.do_diff_percent = 0;
    opts_pls.target_pls_format = 0; % Convert pls to match this format!
    opts_pls.collapse_pls_to_days = 0;
        

    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    
    opts_perm.split_plusminus = 0;
    if groupmode==3 || groupmode==4
        opts_perm.split_plusminus = 1;
    end
    opts_perm.alpha0 = 0.05;
    opts_perm.alpha_bh0 = 0.2;
    %opts_perm.alpha_bh0 = 0.05;
    
    % Map sp parameters
    sp_threshold = 10;
    
    % Group options
    do_group_collapse_pls2days = 0;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dcs.stats_mode = 0;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    opts_PSC.do_stats = 0;
    opts_PSC.do_diagtest = 0;
    opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    opts_PSC.opts_diagtest.threshold = 1.0; 
   

%% Plot switches
plot_on_spect = 0;
plot_on_spectrogram = 0;
plot_on_scatter = 0;
plot_on_bargraph = 0;
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
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm,opts_exclude);
    sp = out_perm.sig_cells;

    % Load bads perm
    [bad_any_perm] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);

    % Map sp's as needed
    [sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any,sp_threshold);

    % Create group template
    mycrit = [2*ones(1,size(sp,2))];
    grt = group0(1);
    grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
    
    switch groupmode
        case 1
            % Run a simple test
            clear group
            i=0;
            i=i+1; group(i)=grt; group(i).criteria(1:2)=[1 2]; group(i).ctgs=1;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:2)=[0 2]; group(i).ctgs=1;   % Ctg1-2 non-deciders
            
            % Calculate legend entries
            group = group.query_legend(group0);
            
        case 2
            % Run a simple test
            clear group
            i=0;
            i=i+1; group(i)=grt; group(i).criteria(1:2)=[2 1]; group(i).ctgs=2;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:2)=[2 0]; group(i).ctgs=2;   % Ctg1-2 non-deciders
            
            % Calculate legend entries
            group = group.query_legend(group0);

        case 3
            % Run a simple test
            clear group
            i=0;
            i=i+1; group(i)=grt; group(i).criteria(1:4)=[1 0 2 2]; group(i).ctgs=1;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 1 2 2]; group(i).ctgs=1;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 0 2 2]; group(i).ctgs=1;   % Ctg1-2 deciders
            
            % Calculate legend entries
            group = group.query_legend(group0);
            
            i=0;
            i=i+1; group(i).legend = 'Sig. Cat';
            i=i+1; group(i).legend = 'Sig. Dog';
            i=i+1; group(i).legend = 'Not Sig.';
            
        case 4
            % Run a simple test
            clear group
            i=0;
            i=i+1; group(i)=grt; group(i).criteria(1:4)=[2 2 1 0]; group(i).ctgs=2;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:4)=[2 2 0 1]; group(i).ctgs=2;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:4)=[2 2 0 0]; group(i).ctgs=2;   % Ctg1-2 deciders

            % Calculate legend entries
            group = group.query_legend(group0);
            
            i=0;
            i=i+1; group(i).legend = 'Sig. Goc';
            i=i+1; group(i).legend = 'Sig. Tad';
            i=i+1; group(i).legend = 'Not Sig.';
            
    end
    
    % Do swap if necessary
    if any(groupmode == 1:2) && swap_mode > 0
        
        switch swap_mode
            case 0;
                swap_map = [1,2,3,4; [1,2,3,4] ];   % Identity
            case 1; swap_map = [1,2; [3,4] ];       % Irrelevant instead of relevant
            case 2; swap_map = [1,2; [1,3] ];       % Ctg1 Rel vs Irrel
            case 3; swap_map = [1,2; [2,4] ];       % Ctg2 Rel vs Irrel
            case 4; swap_map = [1 2; [9 10] ];             % schA schB
            case 5; swap_map = [1:4; 5:8 ];             % FOR CTGSETLI MODE 7
            case {1:5}
                swap_map_pref = unique(round(swap_map/2)','rows')'; 
        end

        if swap_mode == 6
            swap_map = [1,2; 3,4];
            swap_map_pref = swap_map;
        end

        
        group = group_swap.crit_move (group,swap_map_pref);      % THIS BUGGERS UP WHEN USING CTG 11
        group = group_swap.ctgs (group,swap_map);
       
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
    if get_iscromer && groupmode == 0
        if length(group) > 5
            inds=[1:4,9,10];
        else
            inds=[1,2,5];
        end
    else
        inds = 1:N;
%         inds = [1:4];
%         inds = [1:2,5];
    end
            
    
    
    if N <= 5 || ~do_split
        figure;
        [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    else
        %%
        inds = 1:floor(N/2); figure; [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
        inds = floor(N/2)+1:N; figure; [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    end
    
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