
function [wrkspc_buffer, out] = Fg_5_00b_generalized_updated(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)


%% Pull in variables from s

    s = fix_freqband_format(s);     % Converts s.freqband_stats from [18] to [18 18]
    vars_pull(s);
    
    % Annoying stuff that needs to be calculated (parameters based on other parameters...)
    opts_pls.timeband_stats = timeband_stats;
    opts_pls.tf_label_stats = tf_label_stats;
    opts_perm.timeband_perm = timeband_perm;
    opts_perm.tf_label_perm = tf_label_perm;
    is_spectrogram = get_is_spectrogram(sfc_mode) && ~opts_pls.spectrogram2spectra_timeslice;
    
    


%% Import data if running as a function; define dependent parameters

%     running_as_script = ~is_calledby;
%     if ~running_as_script
%         pop_workspace(false);
%     end

        
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
datenum_numeric = out_pls.datenum_numeric;
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
    switch groupmode
        case 0
            % Use default grouping (all pairs, enumerate over ctgs)
            group = group0;
            sp = ~bad_any(:);
            
            % Load legends and colourmap.
            group = group.query_legend(group0);
            if do_custom_colourmap
                if opts_pls.perm2pls || opts_pls.do_diff; group = group.query_colourmap_desired(group0,[0,0,0]);
                else group = group.query_colourmap_desired(group0);
                end
            end

        case 1
            % Load sp's
            [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_perm,bad_any,opts_perm,opts_exclude);
            sp = out_perm.sig_cells;
            
            % Load bads perm
            [bad_any_perm] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);
            
            % Map sp's as needed
            [sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any,sp_threshold);

            % Create group template
            mycrit = [2*ones(1,size(sp,2))];
            grt = group0(1);
            grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
            grt.freqband_perm = freqband_perm; grt.timeband_perm = opts_perm.timeband_perm; grt.tf_label_perm = opts_perm.tf_label_perm;
            grt.Nwind_perm = out_perm.os.Nwind;
            grt.tapers_perm = out_perm.os.params.tapers;
            grt.full_bandwidth_perm = out_perm.os.f_fullbandwidth;
            
            i1 = 1;
            i2 = 2;
            if strcmp(tf_labels_perm(3),'A')        % Select cells based on SchA
                i1 = 1;
                i2 = 1;
            elseif strcmp(tf_labels_perm(3),'B')    % Select cells based on SchB
                i1 = 2;
                i2 = 2;
            end

            % Run a simple test
            clear group
            i=0;
            i=i+1; group(i)=grt; group(i).criteria(i1)=[1]; group(i).ctgs=1;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(i1)=[0]; group(i).ctgs=1;   % Ctg1-2 non-deciders
            i=i+1; group(i)=grt; group(i).criteria(i2)=[1]; group(i).ctgs=2;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(i2)=[0]; group(i).ctgs=2;   % Ctg1-2 non-deciders
%             i=i+1; group(i)=grt; group(i).criteria(5)=[1]; group(i).ctgs=5;   % Ctg1-2 non-deciders
%             i=i+1; group(i)=grt; group(i).criteria(5)=[0]; group(i).ctgs=5;   % Ctg1-2 non-deciders

            if strcmp(tf_labels_perm(3),'A')        % Select cells based on SchA
                group = group(1:2);
            elseif strcmp(tf_labels_perm(3),'B')        % Select cells based on SchB
                group = group(3:4);
            end
            
            % Calculate legend entries
            group = group.query_legend(group0);

        case 2
            % Load sp's
            [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_perm,bad_any,opts_perm,opts_exclude);
            sp = out_perm.sig_cells;

            % Create group template
            mycrit = [2*ones(1,size(sp,2))];
            grt = group0(1);
            grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
            grt.freqband_perm = freqband_perm; grt.timeband_perm = opts_perm.timeband_perm; grt.tf_label_perm = opts_perm.tf_label_perm;
            grt.Nwind_perm = out_perm.os.Nwind;
            grt.tapers_perm = out_perm.os.params.tapers;
            grt.full_bandwidth_perm = out_perm.os.f_fullbandwidth;

            % Run a simple test
            clear group
            i=0;
            i=i+1; group(i)=grt; group(i).criteria(1:2)=[1 0]; group(i).ctgs=1;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:2)=[1 0]; group(i).ctgs=2;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:2)=[0 1]; group(i).ctgs=3;   % Ctg3-4 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:2)=[0 1]; group(i).ctgs=4;   % Ctg3-4 deciders
            
            % Calculate legend entries
            group = group.query_legend(group0);
            
        case 3
            % Load sp's
            [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_perm,bad_any,opts_perm,opts_exclude);
            sp = out_perm.sig_cells;
            
            % Load bads perm
            [bad_any_perm] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);
            
            % Map sp's as needed
            [sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any,sp_threshold);

            % Create group template
            mycrit = [2*ones(1,size(sp,2))];
            grt = group0(1);
            grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
            grt.freqband_perm = freqband_perm; grt.timeband_perm = opts_perm.timeband_perm; grt.tf_label_perm = opts_perm.tf_label_perm;
            grt.Nwind_perm = out_perm.os.Nwind;
            grt.tapers_perm = out_perm.os.params.tapers;
            grt.full_bandwidth_perm = out_perm.os.f_fullbandwidth;
            
            i1 = 1;
            i2 = 2;
            if strcmp(tf_labels_perm(3),'A')        % Select cells based on SchA
                i1 = 1;
                i2 = 1;
            elseif strcmp(tf_labels_perm(3),'B')    % Select cells based on SchB
                i1 = 2;
                i2 = 2;
            end

            % Run a simple test
            clear group
            i=0;
            i=i+1; group(i)=grt; group(i).criteria(i1)=[1]; group(i).ctgs=1;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(i1)=[1]; group(i).ctgs=2;   % Ctg1-2 non-deciders
            i=i+1; group(i)=grt; group(i).criteria(i1)=[0]; group(i).ctgs=1;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(i1)=[0]; group(i).ctgs=2;   % Ctg1-2 non-deciders
            i=i+1; group(i)=grt; group(i).criteria(i2)=[1]; group(i).ctgs=3;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(i2)=[1]; group(i).ctgs=4;   % Ctg1-2 non-deciders
            i=i+1; group(i)=grt; group(i).criteria(i2)=[0]; group(i).ctgs=3;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(i2)=[0]; group(i).ctgs=4;   % Ctg1-2 non-deciders
%             i=i+1; group(i)=grt; group(i).criteria(5)=[1]; group(i).ctgs=5;   % Ctg1-2 non-deciders
%             i=i+1; group(i)=grt; group(i).criteria(5)=[0]; group(i).ctgs=5;   % Ctg1-2 non-deciders

            if strcmp(tf_labels_perm(3),'A')        % Select cells based on SchA
                group = group(1:4);
            elseif strcmp(tf_labels_perm(3),'B')        % Select cells based on SchB
                group = group(5:8);
            end
            
            
            
            
            
            
            
            
            
            % Calculate legend entries
            group = group.query_legend(group0);
            
        case 4
            % Load sp's
            [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_perm,bad_any,opts_perm,opts_exclude);
            sp = out_perm.sig_cells;
            
            % Load bads perm
            [bad_any_perm] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);
            
            % Map sp's as needed
            [sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any,sp_threshold);

            % Create group template
            mycrit = [2*ones(1,size(sp,2))];
            grt = group0(1);
            grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
            grt.freqband_perm = freqband_perm; grt.timeband_perm = opts_perm.timeband_perm; grt.tf_label_perm = opts_perm.tf_label_perm;
            grt.Nwind_perm = out_perm.os.Nwind;
            grt.tapers_perm = out_perm.os.params.tapers;
            grt.full_bandwidth_perm = out_perm.os.f_fullbandwidth;

            % Run a simple test
            clear group
            i=0;
            i=i+1; group(i)=grt; group(i).criteria(1:4)=[1 0 2 2]; group(i).ctgs=1;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 1 2 2]; group(i).ctgs=1;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 0 2 2]; group(i).ctgs=1;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:4)=[2 2 1 0]; group(i).ctgs=2;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:4)=[2 2 0 1]; group(i).ctgs=2;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:4)=[2 2 0 0]; group(i).ctgs=2;   % Ctg1-2 deciders
%             i=i+1; group(i)=grt; group(i).criteria(1:4)=[0,0,0,0]; group(i).ctgs=2;   % Ctg1-2 deciders

%             i=i+1; group(i)=grt; group(i).criteria(3:4)=[1 2]; group(i).ctgs=3;   % Ctg1-2 deciders
%             i=i+1; group(i)=grt; group(i).criteria(3:4)=[1 2]; group(i).ctgs=4;   % Ctg1-2 deciders

%             i=i+1; group(i)=grt; group(i).criteria(3:4)=[2 1]; group(i).ctgs=3;   % Ctg1-2 deciders
%             i=i+1; group(i)=grt; group(i).criteria(3:4)=[2 1]; group(i).ctgs=4;   % Ctg1-2 deciders

%             for i = 1:length(group); group(i).ctgs = 9; end

            % Calculate legend entries
            group = group.query_legend(group0);
            
        case {5,6}      % Different group for each day!
            ids = cellfun(@fname_to_filedate,funames1D);    % ID's
            idsu = unique(ids);                             % ID's unique
            Npairs = length(ids);
            Ndays = length(idsu);
            
            % Calculate SPs for each day
            ids_temp = repmat(ids(:),[1,Ndays]);
            idsu_temp = repmat(idsu(:)',[Npairs,1]);
            sp = ids_temp == idsu_temp;
            
            
            % Create group template
            mycrit = [];
            grt = group0(1);
            grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
            grt.freqband_perm = freqband_perm; grt.timeband_perm = opts_perm.timeband_perm; grt.tf_label_perm = opts_perm.tf_label_perm;
            
            % Set up group criteria
            crit_all = logical(eye(length(idsu)));
            for i = 1:Ndays
                if groupmode == 5
                    group(i)=grt; group(i).criteria = crit_all(i,:); group(i).ctgs = 1; group(i).legend = [funames1D{find(sp(:,i),1,'first')} ' (' num2str(i) ')'];
                elseif groupmode == 6
                    group(i)=grt; group(i).criteria = crit_all(i,:); group(i).ctgs = 2; group(i).legend = [funames1D{find(sp(:,i),1,'first')} ' (' num2str(i) ')'];
                end
            end
            
            % Load legends and colourmap.
            %group = group.query_legend(group0);
            if do_custom_colourmap
                if opts_pls.perm2pls || opts_pls.do_diff; group = group.query_colourmap_desired(group0,[0,0,0]);
                else group = group.query_colourmap_desired(group0);
                end
            end
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
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2,is_spectrogram,opts_pls.timeband_stats);

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
opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars', double(length(group) < 10)};
if plot_on_spect && ~is_spectrogram
    
    disable_split = 0;

    N=length(group);
            
    if (N > 5) && ~disable_split
        %%
        curr_subplots = 1;
        hsp = [];
        
        inds = 1:floor(N/2);
        [hsp, curr_subplots,returns] = new_subplot(2,curr_subplots,hsp,1);    % Create a new subplot entry or start a new figure if current fig is full.
        [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
        
        inds = floor(N/2)+1:N;
        [hsp, curr_subplots,returns] = new_subplot(2,curr_subplots,hsp,1);    % Create a new subplot entry or start a new figure if current fig is full.
        [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
        
    elseif (groupmode == 3) && ~disable_split
        %%
        curr_subplots = 1;
        hsp = [];
        for i = 1:floor(N/2)
            inds = i*2-1:i*2;
            [hsp, curr_subplots,returns] = new_subplot(2,curr_subplots,hsp,1);    % Create a new subplot entry or start a new figure if current fig is full.
            [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
        end
        
    else
        if get_iscromer && groupmode == 0
            if length(group) > 5
                inds=[1:4,9,10];
            else
                inds=[1,2,5];
            end
        else
            inds = 1:N;
        end
        figure;
        [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    end
    
    
%     gr_sp = group_individualize(group(1));
%     opts_PM3Dcs.do_subplots=1;
%     [h1] = plot_matrix3D_custstruct([],gr_sp,opts_PM3D,opts_PM3Dcs);

    
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
if plot_on_spect && is_spectrogram
    %group = group.arrayload('zlims_desired',[0 3]);
    group = group.arrayload('ylims_desired',[0 100]);
    group = group.arrayload('xlims_desired',[ -1.2 1.62]);
%     group = group.arrayload('xlims_desired',[ -.5 1.62]);
    ind = 1:length(group);
    if groupmode == 0; opts_PM3Dsp.show_text_stats = 0; opts_PM3Dsp.show_text_perm = 0; end
    if groupmode == 0; ind = 1:2; end
    hsp = plot_spectrogram_custstruct(group(ind),opts_PM3Dsp);
    %plot_spectrogram_custstruct(group(38:end),opts_PM3Dsp);
    %out = colormap_meld(get_clist(2,'matrix'),[0,0,0],[0,0,0],get_clist(1,'matrix'),64); colormap(out);

    if groupmode == 0
        %% Plot groups
        for i = ind
            hsp.set_gca(i);
            for j = 1:length(tf_avail)
                fullband_T = group(i).Nwind*get_dt;
                fullband_F = group(i).full_bandwidth;
                center_T = tf_avail(j).timeband;
                center_F = tf_avail(j).freqband;
                %mylabel = [num2str(j) ' ' tf_avail(j).label_short];
                mylabel = [num2str(j) '.'];
                ctgs = group(i).ctgs;
                if (ctgs == 1 && strcmp(tf_avail(j).label_short(3),'A')) || ...
                        (ctgs == 2 && strcmp(tf_avail(j).label_short(3),'B'))
                    add_spectrogram_tf_rectangle (center_T,center_F,fullband_T,fullband_F,mylabel,'k', '-', 1);
                end
            end
        end
    end
    
%     %%
%     group = group.arrayload('ylims_desired',[0 40]);
%     group = group.arrayload('xlims_desired',[ -.1 1.62]);
%     gr_sp = group_individualize(group(1));
%     gr_sp(1).zlims_desired = [];
%     plot_spectrogram_custstruct(gr_sp(:),opts_PM3Dsp);

end


%% Test bargraph

if plot_on_bargraph
    fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
end

%% Plot histogram
if plot_on_histogram
%     figure
    %plot_histplot(group,'plot_type','hist_paired');
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

%% Plot real vs imag

if plot_on_real_vs_imag
    % valid_plot_types = {'real_imag_scatter','angle_cdf','angle_pdf','mag_vs_angle_scatter'};
    h = plot_phi_analysis(group(9),'plot_type','real_imag_scatter');
end

%% Plot on group.cells vs days
if plot_on_cells_vs_days
    opts_PDays.x_axis_mode = 3;
    opts_PDays.as_percent = 1;
    opts_PDays.uniform_zaxis_allgroups = 1;
    
    [h1, out] = plot_days_custstruct(group,funames1D,opts_PDays,datenum_numeric);
end


%% Package outputs

% Package data if this is a function
if is_function(mfilename)
    out.wrkspc_buffer = wrkspc_buffer;
    out.group = group;
    out.bad_any = bad_any;
    out.pls = pls;
    out.pls_stats = pls_stats;
    out.abscissa = abscissa;
    out.mypairs = mypairs;
end

end

