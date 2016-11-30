
%% Stage
% clear
% clc
run_setdefaultfig
addpath(genpath('./funcs_plots_preferred/'));

% Data loading parameters
recalc_md = 0;
recalc_cfc = 0;
recalc_psd = 0;
recalc_psd_permute = 0;
recalc_pac = 0;
recalc_pac_permute = 0;
recalc_bootstrap = 0;
reload_data = 0;
if get_iscromer; stage_range = [2,3];
else stage_range = [2,3];
end
sfc_mode=41.6013101;
perm_mode=41.6013111;
cfc_mode=40.60131;  % (amp-amp coupling)
pac_mode=43.601210;
pac_mode_permute = 43.601311;

% perm_mode = sfc_mode+0.000001;   % Make sure these are aligned

group_mode = 1;
gr_submode = 1;             % Which of the ctgs to plot from the bndry, RT modes (i.e. congruent/incongruent; 50 or 60% morphs




% Load preferred units
    % Stage selection

    curr_stage = 3;
    curr_stage_alt = 3;
    curr_stage_sfc = 3;
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
    permdat2pls = 0;            % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    perm2pls = 1;               % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
        perm2pls_dophi = 0;
        sort_pls = 0;           % Sort into preferred and non-preferred
        do_diff = 0;            % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
    

    opts_PM3Dcs.paperfig_mode=0;
    opts_PM3Dcs.do_sgolay=0;
    opts_PSC.paperfig_mode=0;
    opts_PPS.paperfig_mode=0;

    
[~, sfc_subgroups] = decode_sfc_mode(sfc_mode); [fname_suffix_sfc] = build_sfcmode(sfc_mode, sfc_subgroups);
fprintf(['Chosen PSD mode: ' fname_suffix_sfc ' \n']);
[~, cfc_subgroups] = decode_sfc_mode(cfc_mode);
[fname_suffix_cfc] = build_sfcmode(cfc_mode, cfc_subgroups);

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
    plot_PAC = 1;
    plot_PAC_sens = 0;
    compare_boundary_AvB = 0;
    
% Stats
    plot_SFC_spectra_stats = 0;


% Scatter plots
    plot_regress_alpha_vs_beta = 0;

% Permutation test results
    plot_Sch_sens = 0;


%% Import data if running as a function; define dependent parameters

running_as_script = ~is_calledby;
if ~running_as_script
    pop_workspace(false);
end

[~, perm_subgroups] = decode_sfc_mode(perm_mode); [fname_suffix_perm] = build_sfcmode(perm_mode, perm_subgroups);
    
    
%% Groups Definition

if ~exist('group','var')
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
end

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
% fcfc = s.f(1,:);
% clear Cave

% Load PSD data
buff_fieldname_psd = ['psd_' mode2modename(sfc_mode)];
buff_func = @() func_recalc_mode41 (sfc_mode,stage_range,file_range);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname_psd, buff_func,recalc_psd,reload_data);
s = wrkspc_buffer.(buff_fieldname_psd).(stagename(curr_stage_sfc));
S1ave = s.Cave;
f = s.f(1,:);
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

%% PLS variable setup

switch plotmode
    case 2
        ind = find(fcfc >= mean(freqband_corrcoef),1,'first');
        pls = squeeze(Rtr(ind,:,:,:));
        abscissa = fcfc;
        group(1).data_name = 'CFC';
        
        % if pls is small, extend it to the full 11 dimensions
        if size(pls,3) == 2
            temp = ones(size(pls(:,:,end))) * Inf;
            pls = cat(3,repmat(temp,[1,1,10]),pls(:,:,end));
        end
        
    case 3
        ind = find(fcfc >= mean(freqband_corrcoef),1,'first');
        pls = squeeze(Ptr(ind,:,:,:));
        abscissa = fcfc;
        group(1).data_name = 'CFC P';
        
        % if pls is small, extend it to the full 11 dimensions
        if size(pls,3) == 2
            temp = ones(size(pls(:,:,end))) * Inf;
            pls = cat(3,repmat(temp,[1,1,10]),pls(:,:,end));
        end
    case 1
        pls = S1ave;
        group(1).data_name = 'PSD';
        abscissa = f;
        
        
        % Stupid clause to expand pls with fake Infs if it's shorter than desired size (Nctgs = 11)
        pls = pls_pad(pls,sfc_subgroups);
        
        
        
end

for i=1:length(group); if isempty(group(i).xlims_desired); group(i).xlims_desired = [0 90]; end; end

%% Remove bad cells

% Get bad files
bad_any = false(size(bad_clip));
if exclude_clipping; bad_any = bad_any | bad_clip; end
if exclude_60; bad_any = bad_any | bad_60; end

ismonkeyL_temp=~cellfun(@isempty,(strfind(funames1D,'L')));
if excludeL; bad_any = bad_any | ismonkeyL_temp; end
if excludeO; bad_any = bad_any | (~ismonkeyL_temp); end
clear ismonkeyL_temp

%% Swap in perm_mode Cave1 and Cave2 to pls
if permdat2pls
    pls = permCave2pls_swap(wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
        perm2pls_dophi);
    sfc_subgroups = perm_subgroups; fname_suffix_sfc = fname_suffix_perm;
        
    group = group_enum_ctgs(gr2,size(pls),fname_suffix_perm,[0 80]);

end


%% Swap in permutation results to pls
if perm2pls         % Instead of showing raw SFC values, show % successfully passing permutation test
    
    %Differences normalized by null distribution
    pls = perm2pls_swap(wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
        perm2pls_dophi,perm2pls_do_bh,bad_any);
    sfc_subgroups = perm_subgroups; fname_suffix_sfc = fname_suffix_perm;

    %Rework group definition (cheesy)
    group = group_enum_ctgs_diff(gr2,size(pls),fname_suffix_perm,[0 200]);
    
end


%% Sort pls into preferred and non-preferred
if sort_pls
    for i = 1:floor(size(pls,3)/2)
        first = i*2-1;
        second = i*2;
        [pls(:,:,first),pls(:,:,second)] = sort_pref_nonpref(f,pls(:,:,first),pls(:,:,second),freqband_stats);
    end
end

%% Difference between adjacent PLS
if do_diff  % pls(:,:,2)-pls(:,:,1), etc.
    pls = -1*diff(pls,[],3);
    pls = pls(:,:,1:2:end);
    if perm2pls_dophi; pls = get_circ_diff(pls); end
    
    group = group_enum_ctgs_diff(gr2,size(pls),fname_suffix_sfc, [0 80]);
end

%% Calculate sensitivities

S1curr = wrkspc_buffer.(buff_fieldname_psd).(stagename(curr_stage)).Cave;
f2curr = wrkspc_buffer.(buff_fieldname_psd).(stagename(curr_stage)).f; f2curr = f2curr(1,:);

delta_pow = calc_pls_stats(f2curr,S1curr(:,:,end),[1 5],'do_mean_freq',1);
alpha_pow = calc_pls_stats(f2curr,S1curr(:,:,end),[9 13],'do_mean_freq',1);
beta_pow = calc_pls_stats(f2curr,S1curr(:,:,end),[16 20],'do_mean_freq',1);
gamma_pow = calc_pls_stats(f2curr,S1curr(:,:,end),[30 40],'do_mean_freq',1);
all_pow = calc_pls_stats(f2curr,S1curr(:,:,end),[0 100],'do_mean_freq',1);

delta_pow = delta_pow ./ all_pow;
alpha_pow = alpha_pow ./ all_pow;
beta_pow = beta_pow ./ all_pow;
gamma_pow = gamma_pow ./ all_pow;

if ~do_gamma
    sens = [delta_pow; alpha_pow; beta_pow]';
else
    sens = [delta_pow; alpha_pow; beta_pow; gamma_pow]';
end

S1curr_alt = wrkspc_buffer.(buff_fieldname_psd).(stagename(curr_stage_alt)).Cave;
f2curr_alt = wrkspc_buffer.(buff_fieldname_psd).(stagename(curr_stage_alt)).f; f2curr_alt = f2curr_alt(1,:);

delta_pow = calc_pls_stats(f2curr_alt,S1curr_alt(:,:,end),[1 5],'do_mean_freq',1);
alpha_pow = calc_pls_stats(f2curr_alt,S1curr_alt(:,:,end),[9 13],'do_mean_freq',1);
beta_pow = calc_pls_stats(f2curr_alt,S1curr_alt(:,:,end),[16 20],'do_mean_freq',1);
gamma_pow = calc_pls_stats(f2curr,S1curr(:,:,end),[30 40],'do_mean_freq',1);
all_pow = calc_pls_stats(f2curr_alt,S1curr_alt(:,:,end),[0 100],'do_mean_freq',1);

delta_pow = delta_pow ./ all_pow;
alpha_pow = alpha_pow ./ all_pow;
beta_pow = beta_pow ./ all_pow;
gamma_pow = gamma_pow ./ all_pow;

if ~do_gamma
    sens_alt = [delta_pow; alpha_pow; beta_pow]';
else
    sens_alt = [delta_pow; alpha_pow; beta_pow; gamma_pow]';
end

% Cell dates
if ~get_iscromer
    ismonkeyL=~cellfun(@isempty,(strfind(funames1D,'L')));
else
    ismonkeyL=~cellfun(@isempty,(strfind(funames1D,'li'))); % Or something like that...
end
date = cell2mat({md.datenum_numeric});
Nunits = arrayfun(@(s) length(s.lfp_names),md);
date = replicate_inline(date,Nunits); % Make date be 1x549 instead of 1x79 (one entry for each unit)
date = date-min(date);

dateL = date; dateL = dateL - min(date(ismonkeyL)) + 0; dateL(~ismonkeyL) = 0;
dateO = date; dateO = dateO - min(date(~ismonkeyL)) + 0; dateO(ismonkeyL) = 0;
cellnum = [date(:) dateL(:) dateO(:) ismonkeyL(:) ismonkeyL(:)];



%% Derive new spc and spa
switch preferred_mode
    case 1
        % Uses the default spc, extracted above from run_preferred3.m
    case 2
        % Original threshold method
        thresh = 0.15;
        
        spcd = sens(:,:) > thresh;
        spad = sens_alt(:,:) > thresh;
        spc=spcd;
        spa=spad;
    
    case 3
        % Quantile thresholding method
        quant_thresh = 0.7;

        N = size(sens,1);
        sens_curr = sens; temp = quantile(sens_curr(~bad_any,:),quant_thresh); spcd = [sens_curr(:,:) > repmat(temp,N,1)];
        sens_curr = sens_alt; temp = quantile(sens_curr(~bad_any,:),quant_thresh); spad = [sens_curr(:,:) > repmat(temp,N,1)];
        clear temp
        spc = spcd;
        spa = spad;
    
end


%% Parse sensitivity strings used for scatterplots and regression
if sfc_subgroups(2) == 0 
    pls_stats_temp = calc_pls_stats(abscissa,pls,freqband_stats,'do_mean_ctgs',0,'do_mean_freq',0);
    [sens_sc bad_zeros_sc] = parse_sens(senscstr,sens,sens_alt,spc,spa,pls_stats_temp,cellnum);
    [sens_reg bad_zeros_reg]= parse_sens(sensregstr,sens,sens_alt,spc,spa,pls_stats_temp,cellnum);
    clear pls_stats_temp
end



%% Quartile calculations
% Quartile: Splits the last group into SFC and non-SFC cells
q4_sfc = zeros(size(spc,1),1);



%% Calculate group cells
% Extract units sensitive to given schemes
for i = 1:length(group)
    [group(i).cells, group(i).numcells]= get_grouped_cells(group(i),[spc,spa,q4_sfc, ismonkeyL(:)]);
    %fix error here
end

% Remove bad files
spc(bad_any,:) = false;
spa(bad_any,:) = false;
[group] = remove_bads_cells(bad_any,group);                 % Remove bad cells from cell lists
%[bad_any,src,sra,pls,ns,ev,SPG] = remove_bads_data(bad_any,src,sra,pls,ns,ev,SPG);  % Fill dataset values with nans for cells we know are bad!
% group.gettable('numcells');
% group.gettable('right');



%% Groups Loading Setup
% Generate SFC data for plotting and statistics
for i = 1:length(group)
    switch plotmode
        case {1,2,3}
            [group(i).data, group(i).funames] = get_grouped_data(group(i).ctgs,group(i).cells,pls,funames1D);
            [group(i).datastats, group(i).freqband_stats] = calc_pls_stats(abscissa,group(i).data,freqband_stats,'do_mean_ctgs',1);
    end
    group(i).xdata=abscissa;
end


%% pls_stats
% % Calculate pls stat (note - we want this to go AFTER sorting section)
% % We also want to do this after loading data into groups, so we
% % don't duplicate any manipulation of pls (i.e. convert to complex->real).
pls_stats = calc_pls_stats(abscissa,pls(:,:,end),freqband_stats,'do_mean_ctgs',1);          % 11th column for all trials (need to make sure it's all good trials)


if plot_SFC_spectra
    %% Plots Spectra

    figure;
    opts_PM3D = {'do_mean',domean,'do_zscore',do_zscore,'showErrorbars',showError};
    [h1, out.PM3Dcs] = plot_matrix3D_custstruct(abscissa,group,opts_PM3D,opts_PM3Dcs);

    switch plotmode
        %case {2,3}; xlabel('Freq (Hz)'); xlim([0 70]); if do_zscore; ylabel('SFC - normalized'); else ylabel('SFC'); end;
        %case 1; xlabel('Freq (Hz)'); xlim([0 80]); ylabel('PSD');
        
    end
end


if plot_SFC_spectra_stats
    %% Plot Stats
    figure; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
    h
    p
end


if plot_Sch_sens
    [out.PPS] = plot_permute_scatter(...
        wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
        [], ...             % No bootstrap data for this so far
        bad_any,freqband_stats, opts_PPS);
end



if 0
    keyboard
    %% Survey CFC
    
    
    %% Percentage significant
    P = Ptr(:,:,:,end);
    R = Rtr(:,:,:,end);
    alpha = 0.05;
    
    
    
    ind = P < alpha;
    ind = mean(ind,3);
    ind = ind * 100;
    
    figure; imagesc(ind); set(gca,'YDir','normal'); colorbar;
    xlabel('Freq (Hz)'); ylabel('Freq (Hz)')
    
    %% Mean correlation
    figure; imagesc(f,f,mean(R,3)); set(gca,'YDir','normal'); colorbar;
    xlabel('Freq (Hz)'); ylabel('Freq (Hz)')
    
    %% Alpha vs Beta correlations
    alpha_ind = 5;
    beta_ind = 9;
    
    cfcs = squeeze(R(alpha_ind,beta_ind,:));
    cfcs_p = squeeze(P(alpha_ind,beta_ind,:));
    cfcs_h = cfcs_p < 0.05;
    
    figure;
    [n,b] = hist(cfcs,20);
    [n2,b2] = hist(cfcs(cfcs_h),b);
    n = n - n2;                         % Subtract out n2 since it's stacked
    bar([b;b2]',[n;n2]','stacked');
    title(['Correlation between ' num2str(f(alpha_ind)) ' Hz and ' num2str(f(beta_ind)) ' Hz bands']);
    xlabel('Correlation Coefficient');
    ylabel('N cells');
    legend('Non-significant','Significant')
    
    
    %% Plot individuals
%     ind = P < alpha & R < 0;
%     ind2 = find(squeeze(ind(5,9,:)));
%     figure; 
%     for i = 1:length(ind2)
%         imagesc(R(:,:,ind2(i)));set(gca,'YDir','normal'); colorbar;
%         pause
%     end
%     
%     figure; plot(squeeze(R(4,9,ind2)))

end


if 0
    %% Plot PSDs associated with anticorrelated alpha/beta CFC
    
    freq1_ind = 5;  % alpha
    freq2_ind = 9;  % low-beta
    
    P = Ptr(:,:,:,end);
    R = Rtr(:,:,:,end);
    alpha = 0.5;
    
    ind = P < alpha & R < 0;
    ind = squeeze(ind(freq1_ind,freq2_ind,:));
    
    temp = find(f2 < 50);
    figure; plott_ani(f2(temp),S1ave(temp,ind,end));
    
    %figure; plott_ani(f2(temp),S1ave(temp,:,end));
end


if 0
    %% Identify delta, alpha, and beta peaks with findpeaks (not used - powers, below, are better)
    plot_debug = 1;
    plot_on = 1;
    delta_ind = find(f2 >= 1 & f2 <= 6);
    alpha_ind = find(f2 >= 9 & f2 <= 13);
    beta_ind = find(f2 >= 16 & f2 <= 20);
    
    Npsds = size(S1ave,2);
    delta_tr = [];
    alpha_tr = [];
    beta_tr = [];
    for i = 1:Npsds
        ind = findpeaks(S1ave(:,i,end));
        ind = ind.loc;
        
        if plot_debug
            clf; plot(f2,S1ave(:,i,end));
            hold on; plot(f2(ind),S1ave(ind,i,end),'rx')
            xlim([0 40])
            pause
        end
        delta_tr = [delta_tr, ~isempty(find(ind >= delta_ind(1) & ind <= delta_ind(end)))];
        alpha_tr = [alpha_tr, ~isempty(find(ind >= alpha_ind(1) & ind <= alpha_ind(end)))];
        beta_tr = [beta_tr, ~isempty(find(ind >= beta_ind(1) & ind <= beta_ind(end)))];
    end
    
    if plot_on
        figure; plot
    end
end



if 0
    keyboard
    %% Identify delta, alpha, and beta peaks using powers redo
    
    plot_on = 0;
    
    delta_pow = calc_pls_stats(f2,S1ave(:,:,end),[1 5],'do_mean_freq',1);
    alpha_pow = calc_pls_stats(f2,S1ave(:,:,end),[9 13],'do_mean_freq',1);
    beta_pow = calc_pls_stats(f2,S1ave(:,:,end),[16 20],'do_mean_freq',1);
    all_pow = calc_pls_stats(f2,S1ave(:,:,end),[0 100],'do_mean_freq',1);
    
    delta_pow = delta_pow ./ all_pow;
    alpha_pow = alpha_pow ./ all_pow;
    beta_pow = beta_pow ./ all_pow;
    
    figure; plot3(delta_pow,alpha_pow,beta_pow,'.')
    

end


if plot_PAC
    %% PAC Test
        % Test the fraction of electrodes that have significant differences
        % in Phase amplitude coupling
    if floor(pac_mode) == 42
        curr_ctg = 2;       % Ctg 9 - 42 contains i0 = [1,9,11] but 11 is all NaNs due to run_analysis.m error
    elseif floor(pac_mode) == 43
        curr_ctg = 2;       % Ctg 9 - 43 contains i0 = [1,9,10,11] but 11 is all NaNs due to run_analysis.m error
    end
    
    if any(bad_any ~= [bad_60 | bad_clip]); warning('Do not remove any monkeys for this test!'); end
    bad_any = false(1,length(Pac)); warning('Assuming no bad cells')
    Pac_temp = permute(Pac(:,:,curr_ctg,:),[1 4 2 3]);
    pPac_temp = permute(pPac(:,:,curr_ctg,:),[1 4 2 3]);
    sz = size(pPac_temp);
    [sig_cells, p0] = bh_step(reshape(pPac_temp,[sz(1)*sz(2),sz(3)])',0.05,bad_any); 
    sig_cells = reshape(sig_cells,[sz(1),sz(2),sz(3)]);
    
    %figure; imagesc(fdPac1,fdPac2,mean(pdPac_temp(:,:,:) < 0.05,3) > 0.06); colorbar; set(gca,'YDir','normal');  xlabel('Freq Phase (Hz)'); ylabel('Freq Amp (Hz)');
    figure; imagesc(fPac1,fPac2,mean(Pac_temp(:,:,~bad_any),3)); colorbar; set(gca,'YDir','normal');  xlabel('Freq Phase (Hz)'); ylabel('Freq Amp (Hz)'); title('PAC');
    figure; imagesc(fPac1,fPac2,mean(pPac_temp(:,:,~bad_any) < 0.05,3)*100); colorbar; set(gca,'YDir','normal');  xlabel('Freq Phase (Hz)'); ylabel('Freq Amp (Hz)'); title('Percent Elects Significant');
    figure; imagesc(fPac1,fPac2,mean(sig_cells(:,:,~bad_any),3)*100); colorbar; set(gca,'YDir','normal');  xlabel('Freq Phase (Hz)'); ylabel('Freq Amp (Hz)'); title('Percent Elects Significant with 5% FDR');
    
    
end


if plot_PAC_sens
    %% PAC Sensitivity Test
        % Test the fraction of electrodes that have significant differences
        % in Phase amplitude coupling
    if any(bad_any ~= [bad_60 | bad_clip]); warning('Do not remove any monkeys for this test!'); end
    pdPac_temp = permute(pdPac,[1 4 2 3]);
    sz = size(pdPac_temp);
    bad_any = false(1,length(pdPac)); warning('Assuming no bad cells')
    [sig_cells, p0] = bh_step(reshape(pdPac_temp,[sz(1)*sz(2),sz(3)])',0.05,bad_any); 
    sig_cells = reshape(sig_cells,[sz(1),sz(2),sz(3)]);
    
    %figure; imagesc(fdPac1,fdPac2,mean(pdPac_temp(:,:,:) < 0.05,3) > 0.06); colorbar; set(gca,'YDir','normal');  xlabel('Freq Phase (Hz)'); ylabel('Freq Amp (Hz)');
    figure; imagesc(fdPac1,fdPac2,mean(pdPac_temp(:,:,~bad_any) < 0.05,3)*100); colorbar; set(gca,'YDir','normal');  xlabel('Freq Phase (Hz)'); ylabel('Freq Amp (Hz)'); title('Percent Elects Significant');
    figure; imagesc(fdPac1,fdPac2,mean(sig_cells(:,:,~bad_any),3)*100); colorbar; set(gca,'YDir','normal');  xlabel('Freq Phase (Hz)'); ylabel('Freq Amp (Hz)'); title('Percent Elects Significant with 5% FDR');
    temp1 = mean(pdPac_temp(:) < 0.05);
    temp1_ste = std(pdPac_temp(:) < 0.05) / sqrt(length(pdPac_temp(:)));
    dPac9_temp = permute(dPac9,[1 4 2 3]);
    figure; imagesc(fdPac1,fdPac2,mean(dPac9_temp(:,:,:),3)); colorbar; set(gca,'YDir','normal'); xlabel('Freq Phase (Hz)'); ylabel('Freq Amp (Hz)'); title('PAC Ctg1');
    dPac10_temp = permute(dPac10,[1 4 2 3]);
    figure; imagesc(fdPac1,fdPac2,mean(dPac10_temp(:,:,:),3)); colorbar; set(gca,'YDir','normal'); xlabel('Freq Phase (Hz)'); ylabel('Freq Amp (Hz)'); title('PAC Ctg2');
    
    % Plot bar graph comparing average p values for all f1 f2 pairs
%     figure; barwitherr([temp1_ste temp9_ste],[temp1,temp9]); ylim([0.04 0.06]); set(gca,'XTickLabel',{'Ctg1vs2','Ctg9vs10'});ylabel('Mean pValue')

    % Loop through PAC values for all cells
%     for i = 1:length(temp); clf; imagesc(fdPac1,fdPac2,dPac9_temp(:,:,temp(i))); colorbar; set(gca,'YDir','normal');  xlabel('Freq Phase (Hz)'); ylabel('Freq Amp (Hz'); pause; end

end


if 0
    %% Scatter plot of boundary stuff
    pls_stats2 = calc_pls_stats(abscissa,pls,freqband_stats,'do_mean_ctgs',0,'do_mean_freq',0);
    figure; loglog(pls_stats2(:,5),pls_stats2(:,6),'.');
end

if 0
    %% Test relationship between categorization index and differential FFC
    
%     chosen_cells = 
    
    %Sch A Rel
    ind = 9:9+5;
    catFR = statsFR(:,ind); catFR_names = statsnames(ind);
    
    WCD = [catFR(:,1)-catFR(:,2) catFR(:,2)-catFR(:,3) catFR(:,4)-catFR(:,5) catFR(:,5)-catFR(:,6)];
    WCD = abs(WCD);
    WCD = mean(WCD,2);
    
    BCD = mean(abs(catFR(:,3)-catFR(:,4)),2);
    
    CI_SchARel = (BCD-WCD) ./ (BCD+WCD);
    
    
    % Sch B Rel
    ind = 21:21+5;
    catFR = statsFR(:,ind); catFR_names = statsnames(ind);
    
    WCD = [catFR(:,1)-catFR(:,2) catFR(:,2)-catFR(:,3) catFR(:,4)-catFR(:,5) catFR(:,5)-catFR(:,6)];
    WCD = abs(WCD);
    WCD = mean(WCD,2);
    
    BCD = mean(abs(catFR(:,3)-catFR(:,4)),2);
    
    CI_SchBRel = (BCD-WCD) ./ (BCD+WCD);
    
    % Boot up to pairs
    CI_mean = median([CI_SchARel  CI_SchBRel],2);
    
    % Get FFC difference between boundary and non-boundary trials
    pls_stats2 = log(calc_pls_stats(abscissa,pls,freqband_stats,'do_mean_ctgs',0,'do_mean_freq',0));
    pls_boundary_diff = abs(pls_stats2(:,5) - pls_stats2(:,6));
    ulpairs = build_unit_LFP_pairs_all(md)';
    pls_boundary_diff = pls_boundary_diff(ulpairs(:,2));
    
    % Plot FFC difference vs Categorization Index
    figure; plott_fit(CI_mean,pls_boundary_diff,'.');
    
    % Identify a subset of highly responsive electrodes to boundary
    % difference
    ind = quantile(pls_boundary_diff,0.96);
    ind = (pls_boundary_diff > ind);
    
    clear group2
    group2(1:2) = Grp;
    group2(1).cells = ind;
    group2(2).cells = ~ind;
    
    figure; plot_scattergroups(CI_mean,pls_boundary_diff,group2)
    mean(CI_mean(group2(1).cells))
    mean(CI_mean(group2(2).cells))
    
end


if plot_regress_alpha_vs_beta
    %% Test for trends in FFC alpha vs beta peak
        % Found that alpah and beta FFC peaks are positively correlated
    s = wrkspc_buffer.psd_41_60131.(stagename(curr_stage_sfc));
    f_FFC = s.f(1,:);
    Cave_FFC = s.Cave(:,:,end);
    
    [bad_clip, bad_60] = loadall_bad_data_lfp(curr_stage_sfc,file_range,md);
    bad_sfc = bad_clip | bad_60;
    
    ismonkeyL_temp=~cellfun(@isempty,(strfind(s.funames1D,'L')));
    
    if excludeL; bad_sfc = bad_sfc | ismonkeyL_temp; end
    if excludeO; bad_sfc = bad_sfc | ~ismonkeyL_temp; end
    
    bads = bad_sfc;
    
    x = calc_pls_stats(f_FFC,Cave_FFC(:,~bads),[10 12],'do_mean_ctgs',0,'do_mean_freq',0);
    y = calc_pls_stats(f_FFC,Cave_FFC(:,~bads),[16 20],'do_mean_ctgs',0,'do_mean_freq',0);
    xadj = calc_pls_stats(f_FFC,Cave_FFC(:,~bads),[300 500],'do_mean_ctgs',0,'do_mean_freq',0);
    
    figure; plott_fit(x,y,'k.');
    xlabel('FFC Alpha (10-12 Hz)');
    ylabel('FFC Beta (16-20 Hz)');
%     ylabel('FFC Alpha (10-12 Hz)');
    
    %plots_regress([x(:) y(:)],{'SFCbeta','FFCbeta'})
    reg_opts.plot_bargraph=0;
    reg_opts.plot_added=1;
    reg_opts.preproc_transform = 1;
    reg_opts.dependent_var_transform = 1; 
    
    plots_regress([x(:) xadj(:) y(:)],{'Alpha_FFC','adj11','Beta_FFC'},[],[],reg_opts);
    
end

if compare_boundary_AvB
    %% Compare boundary AvB
    
    A_ind = 1;
    B_ind = 2;
    
    opts_PPS_temp = opts_PPS;
    opts_PPS_temp.hide_all_plots = 1;
    opts_PPS_temp.do_bh0=0;
    opts_PPS_temp.sens_chosen_condition = A_ind;
    [out.PPS] = plot_permute_scatter(...
        wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
        [], ...
        bad_any,freqband_stats, opts_PPS_temp);
    schA_pairs = out.PPS.q4_perm;
    
    opts_PPS_temp.sens_chosen_condition = B_ind;
    [out.PPS] = plot_permute_scatter(...
        wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)), ...
        [], ...
        bad_any,freqband_stats, opts_PPS_temp);
    schB_pairs = out.PPS.q4_perm;
    
    if perm2pls
        pls_diff = pls;
    else
        pls_diff = diff(pls,[],3);
        pls_diff = pls_diff(:,:,1:2:end);
    end
    
    mycrit = [schA_pairs(:), schB_pairs(:)];

    
    clear gr_temp;
    gr_template = Grp; gr_template.xlims_desired = [0 200];  gr_template.criteria_alt = []; gr_template.criteria_sfc = []; 
    if perm2pls; gr_template.data_name = 'zscore(\Delta Coh)'; 
    else gr_template.data_name = '\Delta Coh'; 
    end
    i=0;
    
    % Ctg A Sensitive - Pref vs non-preferred
    i=i+1; gr_temp(i) = gr_template; gr_temp(i).criteria=[1 2]; gr_temp(i).ctgs = [A_ind]; gr_temp(i).legend = 'A Es, Ctg A Resp.';
    i=i+1; gr_temp(i) = gr_template; gr_temp(i).criteria=[0 2]; gr_temp(i).ctgs = [A_ind]; gr_temp(i).legend = '~A Es, Ctg A Resp.';
    i=i+1; gr_temp(i) = gr_template; gr_temp(i).criteria=[1 2]; gr_temp(i).ctgs = [B_ind]; gr_temp(i).legend = 'A Es, Ctg B Resp.';
    i=i+1; gr_temp(i) = gr_template; gr_temp(i).criteria=[0 2]; gr_temp(i).ctgs = [B_ind]; gr_temp(i).legend = '~A Es, Ctg B Resp.';
    
    
    % Ctg B Sensitive - Pref vs non-preferred
    i=i+1; gr_temp(i) = gr_template; gr_temp(i).criteria=[2 1]; gr_temp(i).ctgs = [B_ind]; gr_temp(i).legend = 'B Es, Ctg B Resp.';
    i=i+1; gr_temp(i) = gr_template; gr_temp(i).criteria=[2 0]; gr_temp(i).ctgs = [B_ind]; gr_temp(i).legend = '~B Es, Ctg B Resp.';
    i=i+1; gr_temp(i) = gr_template; gr_temp(i).criteria=[2 1]; gr_temp(i).ctgs = [A_ind]; gr_temp(i).legend = 'B Es, Ctg A Resp.';
    i=i+1; gr_temp(i) = gr_template; gr_temp(i).criteria=[2 0]; gr_temp(i).ctgs = [A_ind]; gr_temp(i).legend = '~B Es, Ctg A Resp.';
    
    
    % Both sensitive - AvsB
    i=i+1; gr_temp(i) = gr_template; gr_temp(i).criteria=[1 1]; gr_temp(i).ctgs = [A_ind]; gr_temp(i).legend = 'Multi Es, Ctg A Resp.';
    i=i+1; gr_temp(i) = gr_template; gr_temp(i).criteria=[0 0]; gr_temp(i).ctgs = [A_ind]; gr_temp(i).legend = 'Neither Es, Ctg A Resp.';
    i=i+1; gr_temp(i) = gr_template; gr_temp(i).criteria=[1 1]; gr_temp(i).ctgs = [B_ind]; gr_temp(i).legend = 'Multi Es, Ctg B Resp.';
    i=i+1; gr_temp(i) = gr_template; gr_temp(i).criteria=[0 0]; gr_temp(i).ctgs = [B_ind]; gr_temp(i).legend = 'Neither Es, Ctg B Resp.';
    group2=gr_temp;
    
    
    for i = 1:length(group2); [group2(i).cells, group2(i).numcells]= get_grouped_cells(group2(i),[mycrit]); end
    [group2] = remove_bads_cells(bad_any,group2);                 % Remove bad cells from cell lists
    for i = 1:length(group2)
        [group2(i).data, group2(i).funames] = get_grouped_data(group2(i).ctgs,group2(i).cells,pls_diff,funames1D);
        [group2(i).datastats, group2(i).freqband_stats] = calc_pls_stats(abscissa,group2(i).data,freqband_stats,'do_mean_ctgs',1);
    end
    
    figure; [h1, out.PM3Dcs] = plot_matrix3D_custstruct(abscissa,group2(1:2),opts_PM3D,opts_PM3Dcs);
    figure; [h1, out.PM3Dcs] = plot_matrix3D_custstruct(abscissa,group2(3:4),opts_PM3D,opts_PM3Dcs);
    figure; [h1, out.PM3Dcs] = plot_matrix3D_custstruct(abscissa,group2(5:6),opts_PM3D,opts_PM3Dcs);
    figure; [h1, out.PM3Dcs] = plot_matrix3D_custstruct(abscissa,group2(7:8),opts_PM3D,opts_PM3Dcs);
%     figure; [h1, out.PM3Dcs] = plot_matrix3D_custstruct(abscissa,group2(9:10),opts_PM3D,opts_PM3Dcs);
%     figure; [h1, out.PM3Dcs] = plot_matrix3D_custstruct(abscissa,group2(11:12),opts_PM3D,opts_PM3Dcs);

    
    
end

% Package data if this script is called by a parent function
if is_calledby
    out.wrkspc_buffer = wrkspc_buffer;
    out.group = group;
    out.bad_any = bad_any;
    out.pls = pls;
    out.pls_stats = pls_stats;
end

% Clear everything except wrkspc
clearvars -except wrkspc_buffer out bad_any f f2 pls Cave % Cannot comment this out, because variables will not refresh

