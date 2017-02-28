
function [wrkspc_buffer, out] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,pls_opts)


%% Setup basic parameters
% clc
addpath(genpath('../funcs_plots_preferred/'));

% Data loading parameters
recalc_md = 0;
recalc_sfc = 0;
reload_data = 0;
if get_iscromer; stage_range = [];
else stage_range = [];
end

    if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end
    % SFC mode
    if ~exist('sfc_mode','var'); sfc_mode = 22.4013101; end

    % Stage selection
    if ~exist('curr_stage_sfc','var'); curr_stage_sfc = 2; end
    curr_stage_badclipping = curr_stage_sfc;

    % Unit exclusion
    exclude_nans = opts_exclude.exclude_nans;
    
    % Pls switches
    do_fisher=0;
    if ~exist('freqband_stats','var'); freqband_stats = [16 20]; end
    
    % More pls switches.
    collapse_pls_to_days = pls_opts.collapse_pls_to_days;
    plotmode = pls_opts.plotmode;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    permdat2pls = pls_opts.permdat2pls;             % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    perm2pls = pls_opts.perm2pls;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        perm2pls_do_bh = pls_opts.perm2pls_do_bh;   % Do bh stepup instead of basic test
        perm2pls_dophi = pls_opts.perm2pls_dophi;
        perm2pls_return_mode = pls_opts.perm2pls_return_mode;     % Return mode of perm2pls (4=zscore)
        perm2pls_allow_signed = pls_opts.perm2pls_allow_signed;   % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        perm2pls_split_plusminus = pls_opts.perm2pls_split_plusminus; % 0-Return everything; 1-return pos+negative in separate columns; 2-return only positive (others set to 0); 3-return only negative; 4-return only significant (only makes sense for return_mode=4)
        sort_pls = pls_opts.sort_pls;               % Sort into preferred and non-preferred
        swap_pls = pls_opts.swap_pls;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        do_diff = pls_opts.do_diff;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        do_diff_percent = pls_opts.do_diff_percent;  % Take diff as percent
        do_abs_diff = pls_opts.do_abs_diff;          % Take absolute value after doing diff.
    target_pls_format = pls_opts.target_pls_format;
    timeband_stats = pls_opts.timeband_stats;
    tf_label_stats = pls_opts.tf_label_stats;
    spectrogram2spectra_timeslice = pls_opts.spectrogram2spectra_timeslice;    % If working with a spectrogram, take a slice at time given by timeband_stats.
    spectrogram2ts_freqslice = pls_opts.spectrogram2ts_freqslice;    % If working with a spectrogram, take a slice at time given by freqband_stats.
    spectrogram_normalize_to_baseline = pls_opts.spectrogram_normalize_to_baseline;     % Normalize spectrograms to pre-cue data to a value of 1.0
            spectrogram_baseline_time = pls_opts.spectrogram_baseline_time;             % During pre-cue
            spectrogram_baseline_dolog = pls_opts.spectrogram_baseline_dolog;
    
    
    % Test plot switches
    plot_SFC_spectra = 0;


%% Derive some basic things based on parameters

% Derive suffix name and switches
[~, mode_subgroups] = decode_sfc_mode(sfc_mode);
[fname_suffix, do_adjacent, ctgsetli_mode, thinning_mode, tapers_mode, baseline_subtract, permutation_test, ue_pairs, coh_debias_mode, do_partial_coherence, ctgsetli_mode2, Nwind_mode] = build_sfcmode(sfc_mode, mode_subgroups);
fprintf(['Chosen SFC mode: ' fname_suffix ' \n']);

os = get_options(sfc_mode, curr_stage_sfc);

% Characterize data type;
data_type = mode2datatype(sfc_mode);
is_spectrogram = get_is_spectrogram(sfc_mode);


%% Groups Definition
N_criteria = 1;         % Pick a random number
group = get_group_overrides2(ctgsetli_mode, N_criteria,[],ctgsetli_mode2);
for i = 1:length(group)
    [group(i).xlims_desired, group(i).ylims_desired, group(i).zlims_desired] = get_group_xylims(is_spectrogram && ~spectrogram2spectra_timeslice && ~spectrogram2ts_freqslice, data_type);
    group(i).data_type = data_type;
    group(i).xdata_name = get_group_xdata_name(data_type);
    group(i).freqband_stats = freqband_stats;
    group(i).timeband_stats = timeband_stats;
    group(i).tf_label_stats = tf_label_stats;
    group(i).Nwind = os.Nwind;
    group(i).tapers = os.params.tapers;
    group(i).full_bandwidth = os.f_fullbandwidth;
end
    
%% Load data    


path_matfiles_store = getpath('path_matfiles_store');


% Create buffer variable
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end

% Load Metadata
buff_fieldname = 'currmd';
buff_func = @() func_recalc_md;
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_md,reload_data);
md = wrkspc_buffer.(buff_fieldname).md;


% Estimate file list and file_ranges
% file_range = 1:length(md);
file_range = [];

% Load the SFC data if it's not perm mode
if ~permutation_test
    if strcmp(data_type,'SFC')
        % SFC data
        buff_fieldname = ['sfc_' mode2modename(sfc_mode)];
        buff_func = @() func_recalc_sfc(stage_range,file_range,sfc_mode);
        wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_sfc,reload_data);
    elseif strcmp(data_type,'FFC');
        % Modes 20 and 30 - pairs data - FFC, etc.
        buff_fieldname = ['sfc_' mode2modename(sfc_mode)];
        buff_func = @() func_recalc_pairs(stage_range,file_range,sfc_mode);
        wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_sfc,reload_data);
    elseif strcmp(data_type,'FR');
        % Load FR data
        buff_fieldname = ['sfc_' mode2modename(sfc_mode)];
        buff_func = @() func_recalc_mode41 (sfc_mode,stage_range,file_range);   % We can use this one for FR as well as PSD
        wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_sfc,reload_data);
    else
        % Load PSD data
        buff_fieldname = ['psd_' mode2modename(sfc_mode)];
        buff_func = @() func_recalc_mode41 (sfc_mode,stage_range,file_range);
        wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_sfc,reload_data);
        s = wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc));
    end
end


% Load the Perm data if it is in perm mode
if permutation_test
    perm_mode = sfc_mode;
    buff_mode = perm_mode;
    buff_fieldname = ['per_' mode2modename(buff_mode)];
    buff_func = @() func_recalc_perms(buff_mode,stage_range,file_range);
    wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_sfc,reload_data);
    
end


% Pull in the data
vars_pull(wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sfc)));
% if strcmp(data_type,'FFC');
    vars_pull(wrkspc_buffer.(buff_fieldname).pairsdat.(stagename(curr_stage_sfc)));
% end

if ~exist('phiave','var'); phiave = []; end;



mypairs = mypairs(:,:,end)';
Ncells_shift = Ncells_shift(:);
Ncells_shift = Ncells_shift(:);
Nelects_shift = Nelects_shift(:);
lumap=lumap(:);

% Set up mypairs
if ue_pairs ==4 || ue_pairs == 6    % FFC or PSD
    mypairs = mypairs + [Nelects_shift Nelects_shift];
elseif any(ue_pairs == [0,2,7])     % SFC or FR
    mypairs = mypairs + [Ncells_shift Ncells_shift];
elseif any(ue_pairs == [3])         % Units vs elects (SFC)
    mypairs = mypairs + [Ncells_shift Nelects_shift];
    
end

% For FR time series, shift f so t=0 is aligned with sample on
if strcmp(data_type,'FR')
    temp = get_stagesir(curr_stage_sfc);
    f = f + temp(1)/1e3;
end

% For spectrogram, center f2 time series vector
if is_spectrogram
    %temp = get_stagesir(curr_stage_sfc);
    %f2 = f2 + temp(1)/1e3;
    if exist('stagesir','var')
        if any(f2 > 100)                    % If any time series are greater than 100, we know its in milliseconds.
            f2 = f2/1e3 + stagesir(1)/1e3;  % This works with most recently calculated spectrograms using precalc
        else
            f2 = f2 + stagesir(1)/1e3;      % Annoying case for old spectrogram data
        end
    else
        f2 = f2 - 2.0; warning('Stagesir not found; guessing time range of spectrogram.');
    end
end

%% Remove bad cells

[bad_any,bad_clip,bad_60,fnames_clip] = load_bads(sfc_mode,curr_stage_badclipping,opts_exclude,md,mypairs,funames1D);


%% Import PLS data
if ~permutation_test
    %% Not permutation test simulation
    % Extract pls based on plotmode.
    [pls, group, SPG] = invoke_plotmode(Cave,phiave,group,plotmode,do_fisher,data_type,freqband_stats,f);
    
    
    %% Reshape group & fill in missing data if necessary
        % This is for Cromer data or other weird cases where we don't want
        % to calculate everything pointlessly in analyze_data
    
    if ctgsetli_mode == 0 && size(pls,3) == 4
        group = group([1,9,10,11]);
        pad_mode = 2;
        pls = pls_pad(pls,mode_subgroups,pad_mode);
    elseif ctgsetli_mode == 0 && size(pls,3) == 3
        group = group([1,9,10]);
        pad_mode = 2;
        pls = pls_pad(pls,mode_subgroups,pad_mode);
    elseif ctgsetli_mode == 0 && size(pls,3) == 10
        group = group([1:10]);
        pad_mode = 1;
        pls = pls_pad(pls,mode_subgroups,pad_mode);
    end
    
    


else
    
    % Set up group name
    if ~perm2pls_dophi
        group(1).data_name = ['' get_group_ydata_name(data_type) ''];
    else
        group(1).data_name = ['Angle(' get_group_ydata_name(data_type) ')'];
    end
    
    %% Is permutation test
    if permdat2pls
        %% Swap in perm_mode Cave1 and Cave2 to pls
        
        perm2pls_dophi_temp=0;
        Cave_temp = permCave2pls_swap(wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
            perm2pls_dophi_temp);
        
        perm2pls_dophi_temp=1;
        phiave_temp = permCave2pls_swap(wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
            perm2pls_dophi_temp);
        
        [pls, group, SPG] = invoke_plotmode(Cave_temp,phiave_temp,group,plotmode,do_fisher,data_type,freqband_stats,f);
        clear Cave_temp phiave_temp
        
    end

    if perm2pls         % Instead of showing raw SFC values, show % successfully passing permutation test
        %% Swap in permutation results to pls
        %Differences normalized by null distribution
        pls = perm2pls_swap(wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
            perm2pls_dophi,perm2pls_do_bh,bad_any,perm2pls_return_mode,perm2pls_allow_signed,perm2pls_split_plusminus);
        
        if pls_opts.perm2pls_dophi
            pls(1,:,:) = 1;
            pls(end,:,:) = 1;
        end
        
        if pls_opts.perm2pls_split_plusminus ~= 1
            group = grouppairs_diff(group);
        end


    end
    

end


%% Manipulations to spectrogram data
% Normalize spectrogram to pre-cue data (only works for spectrogram data)
if spectrogram_normalize_to_baseline
    normalize_within_elects = 1;    % 0 - normalize by mean of the pre-cue period.
                                    % 1 - normalize by each electrode individually. This within-electrode
                                    %    normalizing removes electrode-electrode variability.
                                    
                                    
    normalize_across_timerange = 1;
        normalizing_timerange = 0.25;

    pls = unwrap_3Dpls(pls,f,f2);
    
    if ~normalize_across_timerange
        ind = find(f2 >= spectrogram_baseline_time, 1, 'first');     % Find index to normalize along
    else
        ind = find( (f2 >= spectrogram_baseline_time - normalizing_timerange) & (f2 <= spectrogram_baseline_time + normalizing_timerange)); 
    end
    
    if ~normalize_within_elects
        pls = pls ./ repmat(mean(mean(pls(:,~bad_any,:,ind),4),2),[1,size(pls,2),1,size(pls,4)]);
    else
        pls = pls ./ repmat(mean(pls(:,:,:,ind),4),[1,1,1,size(pls,4)]);
    end
    
    if spectrogram_baseline_dolog;
        pls = log(pls);
    end
    
    pls = wrap_3Dpls(pls);
    
end
   
if spectrogram2spectra_timeslice
    if ~is_spectrogram
        error('This options only makes sense if loading spectrogram data.');
    end
    is_spectrogram = 0;     % No longer dealing with spectrograms.
    
    % Take a slice of the spectrogram data
    pls = unwrap_3Dpls(pls,f,f2);
    index = find_closest(f2,mean(timeband_stats));
    pls = pls(:,:,:,index);
    
end

if spectrogram2ts_freqslice
    if ~is_spectrogram
        error('This options only makes sense if loading spectrogram data.');
    end
    is_spectrogram = 0;     % No longer dealing with spectrograms.
    
    % Take a slice of the spectrogram data
    pls = unwrap_3Dpls(pls,f,f2);
    index = find_closest(f,mean(freqband_stats));
    pls = pls(index,:,:,:);
    pls = permute(pls,[4,2,3,1]);
    
    % Frequency now becomes time!
    freqband_stats = timeband_stats;
    f = f2;
    
    for i = 1:length(group);
        group(i).xdata_name = 'Time (s)';
        group(i).xlims_desired = []; group(i).ylims_desired = [];
    end
end

%% If changing target mode
if target_pls_format
    
    % Recursive call to load target mode data.
    pls_opts2 = pls_opts;
    pls_opts2.target_pls_format = 0;
    [wrkspc_buffer, out_target] = load_pls(wrkspc_buffer,target_pls_format,curr_stage_sfc,freqband_stats,opts_exclude,pls_opts2);
    clear pls_opts2
    
    [pls_mapped, bad_any_mapped] = map_pls(sfc_mode,target_pls_format,pls,bad_any,mypairs,out_target.mypairs);
    pls = pls_mapped;
    bad_any = bad_any_mapped;
    funames1D = out_target.funames1D;
    mypairs = out_target.mypairs;
    clear out_target
    
    
    sfc_mode = target_pls_format;
end




%% Shorten groups if needed
if ~exist('pls','var'); error('Pls not found. Likely because were loading perm data, but not in perm2pls or permdat2pls mode.'); end
Nctgs = size(pls,3);
if length(group) > Nctgs
    group = group(1:Nctgs);
end


%% Sort pls into preferred and non-preferred
if sort_pls
    
    for i = 1:floor(size(pls,3)/2)
        first = i*2-1;
        second = i*2;
        [pls(:,:,first),pls(:,:,second)] = sort_pref_nonpref(f,pls(:,:,first),pls(:,:,second),freqband_stats);
    end
end

%% Swap pls entries as desired
if ~isempty(swap_pls)
    source = swap_pls(1,:);
    dest = swap_pls(2,:);
    
    source_pls = pls(:,:,source);
    dest_pls = pls(:,:,dest);
    
    pls(:,:,source) = dest_pls;
    pls(:,:,dest) = source_pls;
    clear source dest source_pls dest_pls
end

%% Difference between adjacent PLS
if do_diff  % pls(:,:,2)-pls(:,:,1), etc.
    if do_diff_percent; pls0 = pls; end
    
    pls = -1*diff(pls,[],3);            % This puts everything in the correct order! (i.e. index 1 - 2)
    pls = pls(:,:,1:2:end);
    
    if do_abs_diff
        pls = abs(pls);
    end
    
    if do_diff_percent
        %warning('This does not work well because dividing by individual spectra can yield singularities. Better to take the mean first and then divide.');
        %keyboard
        
        pls0 = pls0(:,:,1:2:end);
        pls0 = repmat(mean(pls0(:,~bad_any,:,:,:),2),[1,size(pls0,2),1,1,1,1,1]);
        
        pls = pls ./ pls0 * 100;     % Take percent. This divides by non-boundary case for ctgsetli mode4 and 5.
        clear pls0;
    end
    if perm2pls_dophi; pls = get_circ_diff(pls); end

    group = grouppairs_diff(group);     % This is to update the group legends only.
end


%% Remove bad nans

% ID bad data
bad_sfc_nan = get_bad_SFC(pls);

if exclude_nans; bad_any = bad_any | bad_sfc_nan; end


%% Collapse pls to days
if collapse_pls_to_days
    bad_days_thresh = 0.75;
    ismonkeyL_temp = get_ismonkeyL(funames1D);
    sz=size(pls);
    ids = cellfun(@fname_to_filedate,funames1D);    % ID's
    idsu = unique(ids);                             % ID's unique
    Ndays=length(idsu);
    pls_new=zeros([sz(1),Ndays,sz(3)]);
    for i = 1:Ndays
        ind_today=ids==idsu(i);
        pls_new(:,i,:) = mean(pls(:,ind_today & ~bad_any,:),2);
        bad_new(i) = mean(bad_any(ind_today)) > bad_days_thresh;  % More than half of the day's values bad
        bad_60_new(i) = mean(bad_60(ind_today)) > bad_days_thresh;  % More than half of the day's values bad
        bad_clip_new(i) = mean(bad_clip(ind_today)) > bad_days_thresh;  % More than half of the day's values bad
        ismonkeyL_new(i) = mean(ismonkeyL_temp(ind_today)) > bad_days_thresh; if mean(ismonkeyL_temp(ind_today)) < 0.99 && mean(ismonkeyL_temp(ind_today)) > 0.01; warning('expected days monkey identity is not integer. Something wrong'); end
        
    end
    clear ismonkeyL_temp
    bad_sfc_nan = get_bad_SFC(pls_new);
    
    % Reset everythin'
    pls=pls_new;
    bad_any=bad_new;
    bad_60=bad_60_new;
    bad_clip_new=bad_clip;
    ismonkeyL=ismonkeyL_new;
end



if plot_SFC_spectra
%% Plots Spectra

    % Load data into groups
    sp = repmat(~bad_any(:), 1, N_criteria);
    group2 = get_grouped_cells_and_data(group,sp,pls,f,mypairs,bad_any,plotmode,freqband_stats,funames1D);


    figure;
    opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
    inds = 1:length(group2);
    [h1] = plot_matrix3D_custstruct(f,group2(inds),opts_PM3D);
    
end


%% Package output

    % Stupid hacks for spectrogram plots
    if is_spectrogram
        pls = unwrap_3Dpls(pls,f,f2);
    else
        if ~exist('f2','var'); f2 = []; end
    end
    
    % Normalize spectrogram if just doing power
%     if strcmp(group(1).data_name,'PSD')
%         fmat = repmat([1; [f(2:end)]'],[1,size(pls,2),size(pls,3),size(pls,4)]);
%         pls = pls .* fmat.^1.0;
%     end
    
    
    
    out.pls = pls;
    if ~is_spectrogram
        out.pls_stats = calc_pls_stats(f,pls,freqband_stats,'do_mean_ctgs',0);
    else
        out.pls_stats = calc_pls3D_stats(f,f2,pls,freqband_stats,timeband_stats,'ctg_singleton',0);         % Set ctg_singleton to 0 for pls; to 1 for group.data because dimensionality is different!
    end
    out.abscissa = f;
    out.abscissa2 = f2;
    out.bad_any = bad_any;
    out.group = group;
    out.legend = {group.legend};
    out.funames1D = funames1D;
    out.mypairs = mypairs;
    out.sfc_mode = sfc_mode;
    out.datenum_numeric = [md.datenum_numeric];
    out.md = md;

end





function [xlims_desired, ylims_desired,zlims_desired] = get_group_xylims(is_spectrogram, data_type)
    % Sets default Group names
    
    if is_spectrogram                   % Spectrogram
        xlims_desired = [];
        ylims_desired = [];
        zlims_desired = [];
        ylims_desired = [0 100];
        %xlims_desired = [-1.35 2.0];
        %zlims_desired = [0.96 2];
    else
        ylims_desired = [];
        if ~strcmp(data_type,'FR');     % Spectra
            xlims_desired = [0 120];
            zlims_desired = [];
        else                            % Time series (firing rate)
            %xlims_desired = [];
            zlims_desired = [];
            xlims_desired = [-1.3 1.7];
        end
    end

end


function xdata_name = get_group_xdata_name(data_type)

    if ~strcmp(data_type,'FR');
        xdata_name = 'Freq (Hz)';
    else
        xdata_name = 'Time (s)';
    end

end


function ydata_name = get_group_ydata_name(data_type)

    if ~strcmp(data_type,'FR');
        ydata_name = data_type;
    else
        ydata_name = 'Firing Rate';
    end

end

function [pls, group, SPG] = invoke_plotmode(Cave,phiave,group,plotmode,do_fisher,data_type,freqband_stats,f)

        exclude_low_mag_coherence_data = 0;
            lowmag_thresh = 0.4;

        SPG = [];
        switch plotmode
            case {1}
                pls = Cave;
                if do_fisher; pls = fisher_trans(pls); end
                group(1).data_name = get_group_ydata_name(data_type);
                if exist('Ssp','var'); SPG = Ssp; else end
            case 4
                pls = Cave.*exp(1i*phiave);    % Imaginary vector containing amp and phase
                pls = pls_complex2angle_cent0(pls);     % Value between -pi and pi
                pls = abs(pls);                         % Only look at phase or antiphase
                group(1).data_name =  ['Angle(' data_type ')'];
            case 5
                pls = S2ave;
                group(1).data_name =  'SPK PSD';
                SPG = S2AVE;
            case 6
                pls = Cave.*exp(1i*phiave);    % Imaginary vector containing amp and phase
                group(1).data_name =  'ComplexMag';
                fcomplex2real = @abs;
            case 7
                pls = Cave.*exp(1i*phiave);    % Imaginary vector containing amp and phase
                fcomplex2real = @pls_complex2angle_cent0;
                group(1).data_name =  'ComplexPhase';
            case 8
                pls = Cave.*exp(1i*phiave);    % Imaginary vector containing amp and phase
                group(1).data_name =  'FFC';
            case 9
                pls = Cave.*exp(1i*phiave);    % Magnitude of imaginary component of FFC
                pls = abs(imag(pls));
                group(1).data_name =  'iFFC';
        end
        
        
        % Sometimes we may want to exclude data because its coherence is
        % low magnitude. This would mainly be when we're calculating the
        % angle of coherence. (Hard to estimate angle of non-coherent
        % data).
        if exclude_low_mag_coherence_data
            pls_stats = calc_pls_stats(f,Cave,freqband_stats,'do_mean_ctgs',0);
            for i = 1:size(pls_stats,2)
                bads = pls_stats(:,i) < lowmag_thresh;
                pls(:,bads,i) = nan;
            end
        end

end
