
%% Setup parameters
clc
run_setdefaultfig
addpath(genpath('./funcs_plots_preferred/'));
addpath(genpath('./funcs_generalized/'));
clearvars -except wrkspc_buffer popts fv3 myvars*
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end

% function [wrkspc_buffer, out] = Fg_5_generalized_template(wrkspc_buffer,s,opts_exclude,opts_pls,opts_perm)

% Setup params 
    clear group group0

    data_mode = 23;
    switch data_mode
        case 2.0          % SFC
            s.sfc_mode =  2.0511101;
            s.perm_mode = 52.700001001;
        case 2.2          % SFC
            %s.sfc_mode =  2.2511101;
            s.sfc_mode =  2.201711101;
            s.perm_mode = 52.700001001;
        case 2.3          % SFC all pairs
            s.sfc_mode =  2.301711101;
            s.perm_mode = 2.301711101;
        case 3          % SFC spectrogram
            s.sfc_mode =  3.201511101;
            s.perm_mode = 2.201511101;
        case 22         % FFC
            s.sfc_mode =  22.401811101;
            s.perm_mode = s.sfc_mode;
%             s.sfc_mode =  22.401511111;     % Partial FFC
%             s.perm_mode = s.sfc_mode;
%             s.sfc_mode =  22.451411103;     % FFC boundary vs non-boundary
            s.perm_mode = 52.700001001;
        case 23         % FFC spectrogram
            s.sfc_mode =  23.4017101013;
            %s.sfc_mode = 23.4513101031;
            %s.sfc_mode =  23.4018111011;
            s.perm_mode = s.sfc_mode;
            %s.perm_mode = 52.700001001;
            %s.perm_mode =  22.401811101;
        case 41         % PSD
            s.sfc_mode =  41.601811101;
            s.perm_mode = 41.601811101;
            %s.perm_mode = 52.700001001;         % Units
        case 45         % PSD spectrogram
            %s.sfc_mode =  45.6018103043;
            s.sfc_mode =  45.601710101;
            s.perm_mode = s.sfc_mode;
            s.perm_mode = 52.700001001;         % Units
        case 52         % Units time series
            s.sfc_mode  = 52.700001001;
            s.perm_mode = 52.700001001;
        case 52.2         % Units
            %s.sfc_mode  =  52.700201001;
            s.sfc_mode  =  52.700200004;
            s.perm_mode =  52.700001001;
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
    opts_exclude.excludeO = 1; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % Pls switches
    s.freqband_stats = [16 20];
    s.freqband_perm = [16 20];
%     s.freqband_stats = [10 12];
%     s.freqband_perm = [10 12];
%     s.freqband_stats = [20 30];
%     s.freqband_perm = [20 30];
    s.timeband_stats = [.6]; s.timeband_perm = [.6];
    s.tf_label_stats = 'Default'; s.tf_label_perm = 'Default';
    [tf_avail] = get_freqband_timeband(s.perm_mode,opts_exclude); s.tf_avail = tf_avail;
    
%     i=6; s.freqband_stats = tf_avail(i).freqband; s.timeband_stats = tf_avail(i).timeband; s.tf_label_stats = tf_avail(i).label; s.tf_labels_stats = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
%     i=6; s.freqband_perm = tf_avail(i).freqband; s.timeband_perm = tf_avail(i).timeband; s.tf_label_perm = tf_avail(i).label; s.tf_labels_perm = tf_avail(i).label_short; fprintf(['Selecting '  tf_avail(i).label_short ':' tf_avail(i).label '\n']);
    

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
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        opts_pls.do_diff_percent = 0;
        opts_pls.do_abs_diff = 0;            % Take absolute value after doing diff.
    opts_pls.target_pls_format = 0; % Convert pls to match this format!
    opts_pls.collapse_pls_to_days = 0;
    opts_pls.spectrogram2spectra_timeslice = 0;   % If working with a spectrogram, take a slice at time given by timeband_stats.
    opts_pls.spectrogram2ts_freqslice = 0;
    opts_pls.spectrogram_normalize_to_baseline = 0;          % Normalize spectrograms to pre-cue data to a value of 1.0
            opts_pls.spectrogram_baseline_time = -1.199;       % During pre-cue
    
    

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
    s.do_group_normalize_specgram_to_baseline_time = 0;
        s.normalize_within_elects = 1;
        s.specgram_baseline_time = -1.199;
         s.gnsbt_do_log = 0;
    
    % Groupmode
    s.groupmode = 0;   % 0-Use default grouping (all pairs, enumerate over ctgs);
                       % 1:4-Select various subgroups
                       % 5:6-Separate into days
        s.examine_Sch_based_on_animal = 0;          % For animal L, do Cat/Dog; for O do Goc/Tad

    s.swap_mode = 0;

    s.group_do_merge = 1;
        s.groupmerge_operation = 0;

% % Plot switches
    s.plot_on_spect = 1;
    s.plot_on_scatter = 0;
    s.plot_on_bargraph = 0;
    s.plot_on_histogram = 0;
    s.plot_on_real_vs_imag = 0;
    s.plot_on_cells_vs_days = 0;
    s.plot_on_electrode_locations = 0;
    s.plot_on_electrode_distances = 0;
    
    s.doing_cat_vs_dog = (opts_pls.perm2pls && opts_pls.perm2pls_allow_signed == 1) || ...
        (opts_pls.permdat2pls == 1 && opts_pls.do_diff == 1 && opts_pls.do_abs_diff == 0) || ...
        s.group_do_merge == 1;
    
    s.do_custom_colourmap = s.doing_cat_vs_dog;
    s.do_custom_colourmap = false;


% % Plotting options
    paperfig_mode = 1;
    s.opts_PM3Dcs.paperfig_mode=paperfig_mode;
    s.opts_PM3Dcs.stats_mode = 0;
    s.opts_PM3Dcs.do_subplots = 0;
        s.opts_PM3Dcs.max_subplots_per_fig = 16;
    % % Spectrogram plotting options
    s.opts_PM3Dsp.paperfig_mode=paperfig_mode;
    s.opts_PM3Dsp.symmetric_axes = s.doing_cat_vs_dog || opts_pls.spectrogram_normalize_to_baseline || s.do_group_normalize_specgram_to_baseline_time;
    s.opts_PM3Dsp.uniform_zaxis_allgroups = 1;           % Makes z-axis the same for all groups plotted
    s.opts_PM3Dsp.do_subplots = 1;
        s.opts_PM3Dsp.max_subplots_per_fig = 16;
    s.opts_PM3Dsp.show_range_stats = 1;
    s.opts_PM3Dsp.show_range_perm = 1;
        % Overlay Options - Transparency & Contours
        s.PM3Dsp_overlay_opts.do_transparency = 0;
        s.PM3Dsp_overlay_opts.do_contours = 0;
            s.overlay_raw_contours = 0;               % Overlays contours showing raw (non-diffed) FFC values.
            s.swap_in_groupdata_contours = 0;         % Overlays contours showing the same data being plotted in spectrogram (taken from group.data)
            s.swap_in_grouppairs_merge_pvals = 0;     % Overlays contours showing p values
        s.PM3Dsp_overlay_opts.contour_nv = [];
        s.PM3Dsp_overlay_opts.contour_linespec = {'k.'};
        % Stats
        s.PM3Dsp_stats_opts.stats_displaymode = 0;    % 0-no stats; 1 transparency; 2-contours; 3-both (overwrites default overlay settings above)
        s.PM3Dsp_stats_opts.statsfunc = [];
        s.PM3Dsp_stats_opts.stats_comparison = [0];
        s.PM3Dsp_stats_opts.contours_alphas = [0.01 0.001 0.0001];
        s.PM3Dsp_stats_opts.transparency_alpha = [0.01];
    
    s.opts_PSC.paperfig_mode = paperfig_mode;
    s.opts_PSC.remove_dependent = 0;
    s.opts_PSC.do_stats = 1;
    s.opts_PSC.do_diagtest = 0;
    s.opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    s.opts_PSC.opts_diagtest.threshold = 1.0; 
    s.opts_PDays.do_subplots = 1;
        s.opts_PDays.max_subplots_per_fig = 16;
        

%% Pull in variables from s

    s = fix_freqband_format(s);     % Converts s.freqband_stats from [18] to [18 18]
    vars_pull(s);
    
    % Annoying stuff that needs to be updated (parameters based on other parameters...)
    opts_pls.timeband_stats = timeband_stats;
    opts_pls.tf_label_stats = tf_label_stats;
    opts_perm.timeband_perm = timeband_perm;
    opts_perm.tf_label_perm = tf_label_perm;
    is_spectrogram = get_is_spectrogram(sfc_mode) && ~opts_pls.spectrogram2spectra_timeslice && ~opts_pls.spectrogram2ts_freqslice;
    
    [tf_avail] = get_freqband_timeband(s.perm_mode,opts_exclude);


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
md = out_pls.md;
clear out_pls

% Invert pls if needed
[~, sfc_subgroups] = decode_sfc_mode(sfc_mode);
[~, ~, ctgsetli_mode] = build_sfcmode(sfc_mode, sfc_subgroups);
if any(ctgsetli_mode == [4,5]) && opts_pls.do_diff
    pls = pls * -1;
    pls_stats = pls_stats * -1;
end

if overlay_raw_contours
    opts_pls2 = opts_pls;
    opts_pls2.permdat2pls = 1;
    opts_pls2.perm2pls = 0;
    opts_pls2.do_diff = 0;
    [wrkspc_buffer, out_pls2] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls2);
    pls_raw_temp = out_pls2.pls;
    
    % Take average pls_raw between adjacent groups
    sz=size(pls_raw_temp);
    N = floor(sz(3)/2);
    pls_raw = zeros(sz(1),sz(2),N,sz(4));
    for i = 1:N
        pls_raw(:,:,i,:,:) = (pls_raw_temp(:,:,i*2-1,:,:) + pls_raw_temp(:,:,i*2,:,:)) / 2;
    end
    
    PM3Dsp_stats_opts.stats_displaymode = 0;    % Turn off stats!
    PM3Dsp_overlay_opts.do_contours = 1;
    clear opts_pls2 out_pls2 pls_raw_temp
end


%% Build group
% Load sp if necessary
if groupmode ~= 0
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
end

if ~exist('group','var')
    switch groupmode
        case 0                          % All units, all ctgs
            
            group = group0;
            sp = ~bad_any(:);
            
            % Load legends and colourmap.
            group = group.query_legend(group0);
            if do_custom_colourmap
                if opts_pls.perm2pls || opts_pls.do_diff; group = group.query_colourmap_desired(group0,[0,0,0]);
                else group = group.query_colourmap_desired(group0);
                end
            end
        case 1                          % Only sig units, all ctgs
            N = length(group0);
            
            % All ctgs, significant for A
            for i = 1:N
                group(i) = grt;
                group(i).ctgs = i;
                group(i).criteria(1) = [1];
            end
            
            % All ctgs, significant for B
            for i = 1:N
                group(i+N) = grt;
                group(i+N).ctgs = i;
                group(i+N).criteria(2) = [1];
            end
            
            if examine_Sch_based_on_animal
                if opts_exclude.excludeO         % Animal L only
                    group = group(1:N);
                elseif opts_exclude.excludeL     % Animal O only
                    group = group(N+1:2*N);
                end
            end
            
            % Calculate legend entries
            group = group.query_legend(group0);

        case 2                              % Sig vs non-sig
            clear group
            i=0;
            i=i+1; group(i)=grt; group(i).criteria(1)=[1]; group(i).ctgs=1;
            i=i+1; group(i)=grt; group(i).criteria(1)=[0]; group(i).ctgs=1;
            i=i+1; group(i)=grt; group(i).criteria(2)=[1]; group(i).ctgs=2;
            i=i+1; group(i)=grt; group(i).criteria(2)=[0]; group(i).ctgs=2;

            if examine_Sch_based_on_animal
                if opts_exclude.excludeO         % Animal L only
                    group = group(1:2);
                elseif opts_exclude.excludeL     % Animal O only
                    group = group(3:4);
                end
            end
            
            % Calculate legend entries
            group = group.query_legend(group0);
        
        case 3                               % Significant units preferred vs non-preferred
            clear group
            i=0;
            i=i+1; group(i)=grt; group(i).criteria(1)=[1]; group(i).ctgs=1;
            i=i+1; group(i)=grt; group(i).criteria(1)=[1]; group(i).ctgs=2;
            i=i+1; group(i)=grt; group(i).criteria(2)=[1]; group(i).ctgs=2;
            i=i+1; group(i)=grt; group(i).criteria(2)=[1]; group(i).ctgs=1;

            if examine_Sch_based_on_animal
                if opts_exclude.excludeO         % Animal L only
                    group = group(1:2);
                elseif opts_exclude.excludeL     % Animal O only
                    group = group(3:4);
                end
            end
            
            % Calculate legend entries
            group = group.query_legend(group0);


        case 4                                  % For raw data - Ctgs 1-4
            clear group
            i=0;
            i=i+1; group(i)=grt; group(i).criteria(1)=[1]; group(i).ctgs=1;
            i=i+1; group(i)=grt; group(i).criteria(1)=[1]; group(i).ctgs=2;
            i=i+1; group(i)=grt; group(i).criteria(2)=[1]; group(i).ctgs=3;
            i=i+1; group(i)=grt; group(i).criteria(2)=[1]; group(i).ctgs=4;

            if examine_Sch_based_on_animal
                if opts_exclude.excludeO         % Animal L only
                    group = group(1:2);
                elseif opts_exclude.excludeL     % Animal O only
                    group = group(3:4);
                end
            end
            
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
group_temp = group;
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2,is_spectrogram,opts_pls.timeband_stats);

if overlay_raw_contours         % Load raw data too if needed and drop it in data_overlay2 (for contours)
    group_temp = get_grouped_cells_and_data(group_temp,sp,pls_raw,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2,is_spectrogram,opts_pls.timeband_stats);
    for i = 1:length(group)
        group(i).data_overlay2 = group_temp(i).data;
    end
    clear group_temp
    clear pls_raw
end
clear group_temp

if do_group_collapse_pls2days
    group = group_collapse_pls2days(group);
end

if do_group_normalize_specgram_to_baseline_time
    group = group_normalize_specgram_to_baseline_time(group,specgram_baseline_time,normalize_within_elects,gnsbt_do_log);
end


if group_do_merge
    % Flip if ctgsetli mode = 4 or 5.
    N_new=floor(length(group)/2)*2;     % Shorten group if its length is odd
    group = group(1:N_new);
    temp=1:length(group);
    N=length(group);
    %temp=flipud(reshape(1:N,2,floor(N/2)));       % This flips the direction of the diff e.g. making it group2-group1
                                                   % OR (group2-group1) / group1 (for percent diff)
                                                   % Commenting this out makes it group1-group2 (e.g. same direction as pls do_diff)
    temp=temp(:);
    group = grouppairs_merge2(group(temp),opts_pls.perm2pls_dophi,groupmerge_operation);
    
    if swap_in_grouppairs_merge_pvals
        for i = 1:length(group)
            group(i).data_overlay1 = double(group(i).data_pvals < 0.01);
            group(i).data_overlay2 = group(i).data_pvals;
        end
        PM3Dsp_overlay_opts.do_transparency = 1;
        PM3Dsp_overlay_opts.do_contours = 1;
        PM3Dsp_overlay_opts.contour_nv = [0.01,0.001,0.0001];
        PM3Dsp_stats_opts.stats_displaymode = 0;
    end
end

%% Test Plot groups, at last
opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars', double(length(group) < 10)};
if plot_on_spect && ~is_spectrogram
    
    disable_split = 0;

    N=length(group);
            
    if (N > 5) && ~disable_split
        
        curr_subplots = 1;
        hsp = [];
        
        inds = 1:floor(N/2);
        %inds = 1:4;
        [hsp, curr_subplots,returns] = new_subplot(2,curr_subplots,hsp,1);    % Create a new subplot entry or start a new figure if current fig is full.
        [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
        
        inds = floor(N/2)+1:N;
        %inds = 5:8;
        [hsp, curr_subplots,returns] = new_subplot(2,curr_subplots,hsp,1);    % Create a new subplot entry or start a new figure if current fig is full.
        [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
        
    elseif (groupmode == 3) && ~disable_split
        
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

%% Test spectrogram
if plot_on_spect && is_spectrogram
    
    %group = group.arrayload('zlims_desired',[0 3]);
    group = group.arrayload('ylims_desired',[0 100]);
    group = group.arrayload('xlims_desired',[ -1.2 1.62]);
%     group = group.arrayload('xlims_desired',[ -.5 1.62]);

    if any(groupmode == [5,6])
        group = group.arrayload('ylims_desired',[0 60]);
        group = group.arrayload('xlims_desired',[ -1.2 1.62]);
    end
    ind = 1:length(group);
    if groupmode == 0; opts_PM3Dsp.show_text_stats = 1; opts_PM3Dsp.show_text_perm = 1; end
    if any(groupmode == [0,1]); ind = [1,2]; end
    
    
    if swap_in_groupdata_contours
        % Overlays the contours associated with group.data
        PM3Dsp_stats_opts.stats_displaymode = 0;
        PM3Dsp_overlay_opts.do_contours = 1;
        for i = 1:length(group); group(i).data_overlay2 = group(i).data; end
    end
    
    hsp = plot_spectrogram_custstruct(group(ind),opts_PM3Dsp,PM3Dsp_overlay_opts,PM3Dsp_stats_opts);
    %plot_spectrogram_custstruct(group(38:end),opts_PM3Dsp);
    %out = colormap_meld(get_clist(2,'matrix'),[0,0,0],[0,0,0],get_clist(1,'matrix'),64); colormap(out);

    if groupmode == 0
        % Plot groups
        if opts_PM3Dsp.show_range_stats
            subplot_ind = 0;
            for i = ind
                subplot_ind = subplot_ind + 1;
                %hsp.set_gca(subplot_ind);              % For subplotsq
                subplotsq_tight(length(ind),subplot_ind);     % For subplot
                for j = 1:length(tf_avail)
                    fullband_T = group(i).Nwind*get_dt;
                    fullband_F = group(i).full_bandwidth;
                    center_T = tf_avail(j).timeband;
                    center_F = tf_avail(j).freqband;
                    %mylabel = [num2str(j) ' ' tf_avail(j).label_short];
                    mylabel = [num2str(j) '.'];
                    mylegend = group(i).legend;
%                     if (ctgs == 1 && strcmp(tf_avail(j).label_short(3),'A')) || ...
%                             (ctgs == 2 && strcmp(tf_avail(j).label_short(3),'B'))
                    if any(strcmp_substr(mylegend,{'Cat','Dog'})) && any(strcmp_substr(tf_avail(j).label,{'Cat','Dog'})) || ...
                            any(strcmp_substr(mylegend,{'Goc','Tad'})) && any(strcmp_substr(tf_avail(j).label,{'Goc','Tad'})) 
                        add_spectrogram_tf_rectangle (center_T,center_F,fullband_T,fullband_F,mylabel,'k', '-', 1);
                    end
                end
            end
        end
    end
    
    if 0
        %%
        gr_sp = group_individualize(group(1));
        gr_sp = gr_sp.arrayload('ylims_desired',[0 100]);
%         gr_sp = gr_sp.arrayload('xlims_desired',[ -.1 1.62]);
        gr_sp(1).zlims_desired = [];
        plot_spectrogram_custstruct(gr_sp(:),opts_PM3Dsp,PM3Dsp_overlay_opts,PM3Dsp_stats_opts);
    end

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

%% Plot electrode locations
if plot_on_electrode_locations
%     [h1] = plot_electrode_locations(group([1,2]),md,funames1D);
%     [h2] = plot_electrode_locations(group([3,4]),md,funames1D);
    [h1] = plot_electrode_locations(group,md,funames1D);
end

%% Plot electrode distances
if plot_on_electrode_distances
    [hsp] = plot_electrode_distances(group,md,mypairs,bad_any);
    
    compare_distances = 1;
    if compare_distances
        opts_PSC2 = opts_PSC;
        opts_PSC2.remove_dependent = 1;
        fign; [h1] = plot_stats_custstruct(group_out,opts_PSC2);
    end
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

% end

