
function run_analysis (sfc_mode, stage, file_range,i0,j0,outname_embed_i0j0)
%% function run_analysis (sfc_mode, stage, file_range,i0,j0,outname_embed_i0j0)

% % % Test Fries method on Kyle's simulated data. Then test Kyle's method
% and compare
% % % 
format compact
format long g

run_setdefaultfig
if ~isdeployed
    addpath(genpath('./funcs_supporting_local'));
    addpath(genpath('./funcs_run_analysis_precalc'));
end

% if ~exist('sfc_mode','var'); sfc_mode = 22.401311101; end
% if ~exist('sfc_mode','var'); sfc_mode = 52.70121010; end
%if ~exist('sfc_mode','var'); sfc_mode = 2.30151010; end
if ~exist('sfc_mode','var'); sfc_mode = 45.6018101001; end
%if ~exist('sfc_mode','var'); sfc_mode = 40.60000100; end
    % sfc_mode = 1 -> fries
    % sfc_mode = 2 -> chronux SFC
    % sfc_mode = 3 -> coherogram (segmented)
    % sfc_mode = 4 -> SFC analysis trial-by-trial
    % sfc_mode = 5 -> unit analysis # NO LONGER WORKS!
        % sfc_mode = 5.1 -> unit analysis, only this version reports the smoothed mean firing rates instad of firing rates per trial
    % sfc_mode = 6 -> raster_plotting
    % sfc_mode = 7 -> Lepage (not yet implemented)
    % sfc_mode = 8 -> spectrogram
    % sfc_mode = 9 -> evoked responses
    % sfc_mode = 12 -> GBC PL
    % sfc_mode = 13 -> GBC Log
    % sfc_mode = 14 -> SFA PL
    % sfc_mode = 15 -> SFA Log
    % sfc_mode = 16 -> SFA PL_Fast
    
    % 20 <= sfc_mode < 29 - LFP-LFP comparisons (FFC)
    % sfc_mode = 22 -> LFP using Chronux
    
    % 30 <= sfc_mode < 40 - SFC all elects and FFC comparisons
    % sfc_mode = 32 -> SFC (all elects) and FFC using chronux
    
    % sfc_mode = 40 -> Amplitude-Amplitude coupling, trial by trial
    % sfc_mode = 41 -> PSD amplitude
    % sfc_mode = 42 -> Cross-frequency coupling ESC method
    % sfc_mode = 43 -> Cross-frequency coupling MI method
    % sfc_mode = 44 -> Cross-frequency coupling CFC method
    % sfc_mode = 52 -> Used for unit firing rate analysis (the non-parametric kind)
    
    % sfc_mode = *.0   -> unit vs main electrode! - not applicable for unit analysis (sfc_mode=5.*)
    % sfc_mode = *.2   -> unit vs instead adjacent electrode!
    % sfc_mode = *.3   -> unit vs all electrodes
    % sfc_mode = *.4   -> all field-field pairs
    % sfc_mode = *.5   -> all spike-spike pairs
    % sfc_mode = *.6   -> all electrodes
    % sfc_mode = *.7   -> units only
    
    % sfc_mode = *.?X  -> Modify ctgsetli:
    %     X=0, default;
    %     X=1, analyse switch trials;
    %     X=2 reaction times match trials;
    %     X=3 reaction times congruent/incongruent;
    %     X=4 analyse boundary trials;
    %     .... (see ctgsetli)
    % sfc_mode = *.??1  -> Use Mikio's scaled SFC
    % sfc_mode = *.??2  -> Use traditional thinning of SPKs
    % sfc_mode = *.???1 -> Use long bandwidth for multi-taper with asymptotic errors
    % sfc_mode = *.???2 -> Use long bandwidth for multi-taper and Jackknife errors
    % sfc_mode = *.???X -> X={0,3,4,5} - X=0, Use tapers=get_tapers; X=3, tapers=[1,1]; X=4, tapers=[3,5]; X=5, tapers=[5,9]; other options too
    % sfc_mode = *.????1 -> Use baseline subtract
    % sfc_mode = *.?????X -> Permutation mode: X=0, no permutations; X=1, permutation test ctg9/10 (Buschman); X=2, bootstrap test if ctg9/10 same
    % sfc_mode = *.??????1 -> Correct for uneven trial number bias
    % sfc_mode = *.??????X -> Correct for uneven trial number bias
    %                         1-Min trials across pairs of i0s single cell; 2-Min trials across pairs of i0s all cells
    %                         3-Min trials across all i0s single cell; 4-Min trials all i0s all cells
    
if ~exist('stage','var'); stage = 4; end
    % stage = 1 - pre sample on
    % stage = 2 - sample stage
    % stage = 3 - delay stage
    % stage = 4 - sample + delay stage
    % stage = 5 - all
    % stage = 6 - Test1 - match or non-match
    % stage = 7 - TDelay - test delay
    % stage = 8 - Test2 - match (if prior was non-match)
    
% i0 = 17;
% j0 = 1;
if ~exist('i0','var'); i0 = [1:2]; end       % There are a max of 11 categories Roy data; There are fewer for Cromer
if ~exist('j0','var'); j0 = [1:2]; end       % There are a max of 35 units per file for Cromer data; This is about 16 for Roy
%if ~exist('j0','var'); j0 = 8; end       % There are a max of 35 units per file for Cromer data; This is about 16 for Roy
if ~exist('outname_embed_i0j0','var'); outname_embed_i0j0 = 0; end      % Adds a suffix to the file name stating which i0 j0 were swept over


if floor(sfc_mode) == 6
    i0 = 5; 
    outname_embed_i0j0 = 0;
end


% path_metadata = getpath('path_metadata');
path_lfp_sample_delay = getpath('path_lfp_sample_delay');
% path_lfp = getpath('path_lfp');


file_list = get_filelist(path_lfp_sample_delay);
Nfiles = length(file_list);

% Load badtrials
bad_trials = prep_badtrials(stage);

if  ~exist('file_range','var')
    file_range = 1:Nfiles;
%    file_range = [54];
%    file_range = [44];
%    file_range = [70];
%     file_range = [2];
%     file_range = [1:Nfiles];
end

for i = file_range
    %fprintf ('Processing file %d of %d, name %s \n',i,file_range(end),file_list{i});
    
    % General SFC/FFC options
    os = get_options(sfc_mode, stage,i0,j0,[],file_list,i);
    
    calc_stats(file_list{i},sfc_mode,stage,bad_trials,i,i0,j0,outname_embed_i0j0, os);
end

end



function calc_stats(filename_orig, sfc_mode, curr_stage, bad_trials, filenum, i0, j0, outname_embed_i0j0, os)
    %% function calc_stats(filename_orig, sfc_mode, curr_stage, bad_trials, filenum, i0, j0, outname_embed_i0j0)
    if ~exist('filename_orig','var'); filename_orig = 'L091906.mat'; end
    if ~exist('sfc_mode','var'); sfc_mode = 2; end
                % sfc_mode = 1 -> fries
                % sfc_mode = 2 -> chronux
                % sfc_mode = 3 -> Lepage (not yet implemented)
                % sfc_mode = 4 -> mtcoherencysegc (segmented)

    % Note: indices refers to indices of lfp matrix
    % Trials refers to indices of the ctgX_trials arrays
    
%     profile on
    
    ticID2 = tic;
    
    % Generate outpath
    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix, do_adjacent, ctgsetli_mode, thinning_mode, tapers_mode, baseline_subtract, permutation_test, ue_pairs, coh_debias_mode, do_partial_coherence, ctgsetli_mode2, Nwind_mode] = build_sfcmode(sfc_mode, mode_subgroups);
    sfc_mode_group = floor(sfc_mode); 
    fprintf(['Chosen SFC mode: ' fname_suffix ' \n']);
    outpath = fullfile(getpath('path_buffer_curr'),['test2mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(curr_stage)]);
    if outname_embed_i0j0
        outpath = fullfile(outpath,'i0j0_split'); [~,filename,~] = fileparts(filename_orig);
        filename = [filename '_i' num2str(i0(1),'%2.2d') '_j' num2str(j0(1),'%2.3d') '.mat'];
    else
        filename = filename_orig;
    end
    mkdir_dave(outpath);
    outname = fullfile(outpath,filename);
    
    if exist(outname,'file')
        fprintf('File %s exists...skipping\n',outname);
        return
    end
    
    % Assign function handle for spike-field or field-field coherence
    coherence_func = @calc_sfc_struct;
    
    % Load data
    load (fullfile(getpath('path_lfp_sample_delay'),filename_orig));
    lfp_sample = double(lfp_sample);
    lfp_sample_indices = double(lfp_sample_indices); % Necessary for parfor, so knows data is loaded
    
    % Load ctgsetli for all files if necessary
    if coh_debias_mode == 2 || coh_debias_mode == 4           % Calculate Ctgsetli for all cells. Only do this if necessary, b/c it's slow.
        tic
        [ctgsetli_A, ctgsetli_good_A] = calc_actual_good_ctgsetli(sfc_mode,bad_trials);
        % save('temp.mat','ctgsetli_A','ctgsetli_good_A');
        et=toc;
        % fprintf(['Elapsed time to calculate all good ctgsetli = ' num2str(et) '.']);
        % load('temp.mat')
        % warning('Uncomment our original calculation.');
        
        
        ctgsall.ctgsetli_A = ctgsetli_A;
        ctgsall.ctgsetli_good_A = ctgsetli_good_A;
        ctgsall.ctgsetli_mins_allcells = calc_ctgsetli_mins_all(ctgsetli_good_A);
    else
        ctgsall=[];
    end
    
    % Do analysis of coherence! (permutation test or otherwise)
    if permutation_test == 0    % Default mode
        permutation_mode = 0;
        sfc = coherence_func(lfp_sample, lfp_sample_indices, filename_orig, sfc_mode, curr_stage, bad_trials, ctgsall, filenum, i0, j0, permutation_mode, os);
        
    elseif (permutation_test == 1) || (permutation_test == 2)                    % Do permutation test for significance
        
        
        if mod(length(i0),2) ~= 0
            error('The permutation test requires an length(i0) to be even. Pairs of ctgs must be compared.');
        end
        
        Ni0pairs = floor(length(i0)/2);
        for pi = 1:Ni0pairs     % i0 pair index
            ticID = tic;
            i0curr = i0(pi*2-1:pi*2);
            
            % Get actual statistic
            permutation_mode = 0;
            sfc_stat = coherence_func(lfp_sample, lfp_sample_indices, filename_orig, sfc_mode, curr_stage, bad_trials, ctgsall, filenum, i0curr(1:2), j0, permutation_mode, os);
            isphi = isfield(sfc_stat{end},'phiave');
            if length(sfc_stat) < i0curr(2); continue; end; % If i0 was larger than size(ctgsetli0), continue
            Cave1 = sfc_stat{i0curr(1)}.Cave;
            Cave2 = sfc_stat{i0curr(2)}.Cave;
            if isphi
                phi1 = sfc_stat{i0curr(1)}.phiave;
                phi2 = sfc_stat{i0curr(2)}.phiave;
            end
            dCave = Cave1 - Cave2;
            %clear sfc_stat

            % Generate null distribution
            sz = size(dCave);
            Npermut = 50;
                if floor(sfc_mode) == 44
                    Npermute = 50;          % Do fewer permutations when in cfc mode for calculating PAC
                end
            dCave_null = zeros([sz, Npermut]);
            Cavep1 = zeros([sz, Npermut]);
            Cavep2 = zeros([sz, Npermut]);
            phip1 = zeros([sz, Npermut]);
            phip2 = zeros([sz, Npermut]);
            permutation_mode = permutation_test;    % Turn on the appropriate permutation test (there is a handy 1:1 map between the test and the mode)
            clear sfc_temp_all sfc_temp
            for k = 1:Npermut
                k
                sfc_temp_all{k} = coherence_func(lfp_sample, lfp_sample_indices, filename_orig, sfc_mode, curr_stage, bad_trials, ctgsall, filenum, i0curr(1:2), j0, permutation_mode, os);
            end
            
            for k = 1:Npermut
                sfc_temp = sfc_temp_all{k};
                Cavep1(:,:,k) = sfc_temp{i0curr(1)}.Cave;
                Cavep2(:,:,k) = sfc_temp{i0curr(2)}.Cave;
                if isphi
                    phip1(:,:,k) = sfc_temp{i0curr(1)}.phiave;
                    phip2(:,:,k) = sfc_temp{i0curr(2)}.phiave;
                end
                %dCave_null(:,:,k) = sfc1p(k).Cave - sfc2p(k).Cave;
            end

            dCave_null = Cavep1 - Cavep2;

            sz = size(dCave_null);
            %pvals = sum(dCave_null >= abs(repmat(dCave,[1,1,sz(3)])),3) ./ (sz(3)*ones([sz(1:2)]));
            %critval = prctile(dCave_null,95,3);

            % Data fields
            %sfc.dCave_null = dCave_null;
            %sfc{1}.dCave = dCave;
            pi2 = floor(i0curr(2)/2);
            sfc{pi2}.Cave1 = Cave1;
            sfc{pi2}.Cave2 = Cave2;
            if isphi
                sfc{pi2}.phi1 = phi1;
                sfc{pi2}.phi2 = phi2;
            end
            sfc{pi2}.Cavep1 = Cavep1;
            sfc{pi2}.Cavep2 = Cavep2;
            if isphi
                sfc{pi2}.phip1 = phip1;
                sfc{pi2}.phip2 = phip2;
            end
            sfc{pi2}.elapsedTime = toc(ticID);
            sfc{pi2}.Ntraces = sfc_temp{i0curr(1)}.Ntraces;
            sfc{pi2}.Ntraces2 = sfc_temp{i0curr(2)}.Ntraces;
            sfc{pi2}.sum_ctgsetli_good1 = sfc_temp{i0curr(1)}.sum_ctgsetli_good;
            sfc{pi2}.sum_ctgsetli_good2 = sfc_temp{i0curr(2)}.sum_ctgsetli_good;            
            sfc{pi2}.sum_ctgsetli_maxes1 = sfc_temp{i0curr(1)}.sum_ctgsetli_maxes;
            sfc{pi2}.sum_ctgsetli_maxes2 = sfc_temp{i0curr(2)}.sum_ctgsetli_maxes;
            
            sfc{pi2}.mypairs = sfc_temp{i0curr(1)}.mypairs;
        end
        
        % Singleton fields
        singletonfieldnames = sfc_temp{1}.singletonfieldnames;
        for i = 1:length(singletonfieldnames)
            if isfield(sfc_temp{1},singletonfieldnames{i})
                sfc{1}.(singletonfieldnames{i}) = sfc_temp{1}.(singletonfieldnames{i});
            end
        end
        
        if ~outname_embed_i0j0
            % Merge permutation tests!
            sfc = permute_null_to_pvals(sfc);
        end
        
    end
    
    toc(ticID2);
%     profile viewer
%     profile off
    save(outname,'sfc');
   
end


function sfc = calc_sfc_struct(lfp_sample, lfp_sample_indices, filename_orig, sfc_mode, curr_stage, bad_trials, ctgsall, filenum, i0, j0, permutation_mode, os)
    %% sfc = calc_sfc_struct(lfp_sample, lfp_sample_indices, filename_orig, sfc_mode, curr_stage, bad_trials, ctgsall, filenum, i0, j0, permutation_mode)
    
    
    % Pull options
    % (alternatively, use vars_pull(os))
    mode_subgroups = os.mode_subgroups;
    sfc_mode_group = os.sfc_mode_group;

    exclude_clipping_trials = os.exclude_clipping_trials;
    drop_clipping_electrodes_for_partial_coh = 1;
    curr_stage = os.curr_stage;
    Nwind = os.Nwind;
    dt = os.dt;
    thresh = os.thresh;
    params = os.params;
    fract_overlap = os.fract_overlap;

    clear os
    
    % Function specific options
    plot_on = 0;
    plot_debug = 0;
    plot_debug1_thinning = 0;
    StarRate = get_rateStar;
    loadall_lfps = 0;
    
    
    % Pull filename info (note - these variables are already pulled in get_options; I am just loading them again for simplicity
    [fname_suffix, do_adjacent, ctgsetli_mode, thinning_mode, tapers_mode, baseline_subtract, permutation_test, ue_pairs, coh_debias_mode, do_partial_coherence, ctgsetli_mode2, Nwind_mode] = build_sfcmode(sfc_mode, mode_subgroups);
    
    if sfc_mode_group == 31 || sfc_mode_group == 52;    % Can't do baseline subtract, since we're using average signals here.
        baseline_subtract = 0;
        warning('Baseline subtract was turned off automatically.\n');
    end
    
    % If doing partial coherence, make sure we're in mode 22
    if do_partial_coherence == 1
        if sfc_mode_group ~= 22
            warning('For partial coherence mode, must be in sfc_mode == 22.*. Wont work with other modes because matrix must be symmetric to take inverse.');
        end
    end

    md = load (fullfile(getpath('path_metadata'),filename_orig));     % Load metadata
    %load (fullfile(getpath('path_lfp'),filename));     % Load full file
    %lfpmat = double(lfpmat_int16); clear lfpmat_int16;
    

    
    if sfc_mode_group ~= 5
        %lfp_indices = recall_lfpir(lfp_sample_indices); % I don't think I use this
        lfp_indices_abs = bounds_to_indices([lfp_sample_indices(:,1) lfp_sample_indices(:,3)]');
        %lfp_indices_abs2 = recall_lfpia(lfp_sample_indices);
    end
    
    Nsamples = length(md.sample_on_times); %Number of "good" samples - when trial reached sample_on stage
    Nunits = length(md.unit_names);    % Number of units tracked
    Nelects = length(md.lfp_names);      % Number of electrodes
    
    
    
    % Set up mypairs for iteration
    switch ue_pairs
        case {0,2};                                             % Current unit versus current/adjacent field (SFC)
            Ncoherences = Nunits;                           
            mypairs = [1:(Nunits); 1:(Nunits)]';
        case 3                                                  % Current unit all fields (SFC)
            % Generate all possible unit-electrode pairs
            [units1,lfp2] = meshgrid(1:Nunits,1:Nelects);
            mypairs = [units1(:) lfp2(:)];
            %mypairs = mypairs(units1 ~= 0,:);
            clear units1 lfp2
            Ncoherences = size(mypairs,1);
        case 4                                                  % Current field versus all field (FFC)
            mypairs = upper_triangle_pairs(Nelects);
            Ncoherences = size(mypairs,1);
        case 5                                                  % Current unit versus all units (SSC)
            mypairs = upper_triangle_pairs(Nunits);
            Ncoherences = size(mypairs,1);
        case 6                                                  % Current electrode only
            mypairs = [1:(Nelects); 1:(Nelects)]';
            Ncoherences = size(mypairs,1);
        case 7                                                  % Current unit only (non-parametric unit analysis)
            mypairs = [1:(Nunits); 1:(Nunits)]';
            Ncoherences = size(mypairs,1);
    end
    
    
    % If we're doing partial coherence, we want to remove the bad
    % electrodes from the pool since they might remove signal components
    % we would otherwise want to keep.
    if do_partial_coherence
        
            % Identify bad clipping electrodes
            path_badchannels = fullfile(getpath('path_badchannels'));
            name_clip = ['id_clipping_s' num2str(curr_stage) '.mat'];
            s_clip = load(fullfile(path_badchannels,name_clip));

            % Optional - use custom thresholds for selecting bad files
            clth = get_clipping_threshold_trials;
            li_clip_all = cat(2,s_clip.clippingstat{:});
            bad_clip_all = li_clip_all > clth(2);   % All files; used for testing purposes only

            li_clip = s_clip.clippingstat{filenum};
            bad_clip = li_clip > clth(2);           % Current file
            
            % If all electrodes are "bad", pick 2 at random and make them
            % good so we don't crash. These will be excluded later.
            while sum(~bad_clip) < 2; bad_clip(unidrnd(length(bad_clip),1,1)) = false; end
            
            bad_clip_pairs = bad_clip(mypairs(:,1)) | bad_clip(mypairs(:,2));
            clear bad_clip_all
            
            
            
        if ~drop_clipping_electrodes_for_partial_coh
            bad_clip = false(1,Nelects);
            bad_clip_pairs = false(1,size(mypairs,1));
        end
    end
            
% % % % % % % % % % % % % % % % % % % % % %     Calculate Ctgsetli's % % % % % % % % % % % % % % % % % % % % % %
    
    % Precalculate ctgsetli's
    [ctgsetli0,N_ctgs_base,N_ctgs_extras] = get_ctgsetli(sfc_mode,md,ctgsetli_mode,ctgsetli_mode2,filenum);
    N_categories = size(ctgsetli0,2);    % Update N_categories
    
    if strcmp(filename_orig,'L112106.mat') && (ue_pairs ~= 5) && (ue_pairs ~= 7)
        ctgsetli0 = ctgsetli0(1:size(lfp_sample,2),:);       % Shorten ctgsetli if we're missing lfp data. Should only be for file #29. Only do this if we are consider LFP
    end
    
    
    if do_partial_coherence
        j0_precalc = 1:Ncoherences;     % Precalc all ctgsetli_goods
    else
        j0_precalc = j0;                % Just precalc the ones of interest
    end
    
    if ~do_partial_coherence
        
        [ctgsetli_good, ctgsetli_maxes] = calc_ctgsetli_goods(j0_precalc,Ncoherences,sfc_mode_group,exclude_clipping_trials,bad_trials,filenum,md,ctgsetli0,ue_pairs,permutation_mode,mypairs,ctgsall,coh_debias_mode,i0);
        
        
    else
    
        % Generate a ctgsetli good set across all electrodes, since we're using
        % all of them.
    
        if ~all(bad_clip_pairs)
            bad_trials_pairs = zeros(size(bad_trials{filenum},1),size(mypairs,1));
            for j = j0_precalc
                if j > Ncoherences; continue; end
                curr_lfp = mypairs(j,1);
                curr_lfp2 = mypairs(j,2);
                bad_trials_pairs(:,j) = bad_trials{filenum}(:,curr_lfp) | bad_trials{filenum}(:,curr_lfp2);
            end
            
            bad_trials_any = any(bad_trials_pairs(:,~bad_clip_pairs),2);
            ctgsetli_good_all = ctgsetli0 & repmat(~bad_trials_any(:),1,size(ctgsetli0,2));
            
            if permutation_mode == 1
                if length(i0) ~= 2; warning('When permutation_mode=1, length(i0) must be 2.'); return; end
                ctgsetli_good_all = permute_ctgsetli(ctgsetli_good_all, i0);
            end
            ctgsetli_maxes_all = calc_ctgsetli_maxes(ctgsetli_good_all,ctgsall,coh_debias_mode);
            
        end
    end
    
    
    
    
% % % % % % % % % % % % % % % % % % % % % %   End  Calculate Ctgsetli's % % % % % % % % % % % % % % % % % % % % % %
    
    
    
    
    % Load data appropriate to current stage
    if any(ue_pairs == [3,4,6])
        [currlfp_all0,lfpia_stage] = get_LFP_ts(lfp_sample,curr_stage,md,lfp_sample_indices);
    end
    if ue_pairs == 3 || ue_pairs == 5 || (ue_pairs == 2 && do_partial_coherence)      % Only calculate this if doing spikes, since profiler says it's slow.
        [currspike_all0] = get_spike_ts_all(Nunits,curr_stage,md);
    end
    
    % Build up unit names
    clear unit_names
    switch ue_pairs
        case {0,2,7}; unit_names = md.unit_names;
        case 3; for kk = 1:Ncoherences; unit_names{kk} = [md.unit_names{mypairs(kk,1)} ':' md.lfp_names{mypairs(kk,2)}]; end
        case 4; for kk = 1:Ncoherences; unit_names{kk} = [md.lfp_names{mypairs(kk,1)} ':' md.lfp_names{mypairs(kk,2)}]; end
        case 5; for kk = 1:Ncoherences; unit_names{kk} = [md.unit_names{mypairs(kk,1)} ':' md.unit_names{mypairs(kk,2)}]; end
        case 6; unit_names = md.lfp_names;
    end

    
    for i = i0
        i;
        
        
        if i > N_categories; continue; end
        switch sfc_mode_group
            case 1
                sfc{i}.Cave = [];       % Make it a cell array to allow for variation amongst # of categories
                sfc{i}.lfp_wind = [];       % Would be useful if we ever want to store individual windows
                sfc{i}.S12ave = [];
                sfc{i}.S1ave = [];
                sfc{i}.S2ave = [];
                sfc{i}.phiave = [];

            case {2, 22, 23, 32}
                sfc{i}.Cave = [];
                sfc{i}.S12ave = [];
                sfc{i}.S1ave = [];
                sfc{i}.S2ave = [];
                sfc{i}.phiave = [];
                sfc{i}.confC = [];
                sfc{i}.phistd = [];
                sfc{i}.Cerr1 = [];
                sfc{i}.Cerr2 = [];

            case 3
                sfc{i}.Ssp = [];
                sfc{i}.S12AVE = [];
                sfc{i}.S1AVE = [];
                sfc{i}.S2AVE = [];
                sfc{i}.ConfC = [];
                sfc{i}.Phistd = [];
                
            case 4
                sfc{i}.S12ave = [];
                sfc{i}.S1ave = [];
                sfc{i}.S2ave = [];

            case 5
                sfc{i}.spikerates = [];
                
            case 8
                sfc{i}.Ssp = [];
            case 9
                sfc{i}.E = [];
            case 12
                sfc{i}.Cave = [];
                sfc{i}.S12ave = [];
                sfc{i}.S1ave = [];
                sfc{i}.S2ave = [];
                sfc{i}.phiave = [];
                sfc{i}.alpha = [];
                sfc{i}.p = [];
            case 15
                sfc{i}.Cave = [];
                sfc{i}.S12ave = [];
                sfc{i}.S1ave = [];
                sfc{i}.S2ave = [];
                sfc{i}.phiave = [];
                sfc{i}.alpha = [];
                sfc{i}.p = [];
                
            case 40
                sfc{i}.Cave = [];   % R value
                sfc{i}.Cerr1 = [];
                sfc{i}.Cerr2 = [];
                sfc{i}.Ptr = [];    % Probability of R value
                
            case 41
                sfc{i}.Cave = [];
                sfc{i}.Cerr1 = [];
                sfc{i}.Cerr2 = [];
                
            case {42,43,44}
                sfc{i}.Cave = [];
                sfc{i}.p = [];
                
            case 52
                sfc{i}.Cave = [];
        end
        
        if sfc_mode_group ~=5
            sfc{i}.Ntraces = [];
            sfc{i}.Ntraces_below_thresh = [];
            sfc{i}.spikerate_mu = [];
            sfc{i}.Ntraces_below_thresh2 = [];
            sfc{i}.spikerate_mu2 = [];
            sfc{i}.adj_was_missing = [];
            sfc{i}.elapsedTime = [];
            sfc{i}.mypairs = [];
            %sfc{i}.sum_ctgsetli2 = [];     % Actually this isn't necessary; we have Ntraces to save this info
            sfc{i}.sum_ctgsetli_good = [];
            sfc{i}.sum_ctgsetli_maxes = [];
        end
        
        
        if do_partial_coherence
            % Get current trials (note this is independent of j now)
            curr_trials = ctgsetli_good_all(:,i);
            Ntrials_min = ctgsetli_maxes_all(i);
            
            % Do baseline subtraction from all electrodes
            if any(ue_pairs == [0,2,3,5])
                currspike_all_ctg = do_baseline_subtract(currspike_all0(:,curr_trials,:),baseline_subtract);
            end
            if any(ue_pairs == [3,4,6])
                currlfp_all_ctg = do_baseline_subtract(currlfp_all0(:,curr_trials,:),baseline_subtract);
            end
            
            if permutation_mode == 2
                warning('Need to test that bootstrap works with partial coherence; I have not tested this yet');
                Ntrials_curr = sum(curr_trials);
                perm = floor(unifrnd(1,Ntrials_curr+1-0.0001,1,Ntrials_curr)); % Uniform sampling from 1 to Ntrials
                if any(ue_pairs == [2,3,5]); currspike_all0 = currspike_all0(:,perm,:); end
                if any(ue_pairs == [3,4,6]); currlfp_all0 = currlfp_all0(:,perm,:); end
            end
            
            % Run partial coherence, precalculating everything in advance!
            params.Fs = 1/dt;
            params.trialave = 1;        
            params.fname = @coherencyc;
            
            % Compute partial coherence only for "good" electrode pairs
            [Cave_all_temp,phiall_temp,~,~,~,f]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@partcoh_stupidwrapper_full,currlfp_all_ctg(:,:,~bad_clip),currlfp_all_ctg(:,:,~bad_clip),params);
            
            % Fill in remaining pairs with NaN's; we should now be back up
            % to full dimensionality
            Cave_all0 = Inf*ones(size(Cave_all_temp,1),Nelects,Nelects);
            Cave_all0(:,~bad_clip,~bad_clip) = Cave_all_temp;
            phiall = Inf*ones(size(phiall_temp,1),Nelects,Nelects);
            phiall(:,~bad_clip,~bad_clip) = phiall_temp;
            clear Cave_all_temp phiall_temp
            
            % Rearrange into 1D based on mypairs.
            sz=size(Cave_all0);
            ind = sub2ind([sz(2),sz(3)],mypairs(:,1),mypairs(:,2));
            Cave_all = Cave_all0(:,ind);
            phiall = phiall(:,ind);
%             Cave_all = Cave_all
        end
        
        
        for j = j0
            j;
            if j > Ncoherences; continue; end
            
            ticID = tic;
            
            
            [i j];
            curr_unit = j;
            unit_names2= strrep(unit_names,'_',' ');
            adj_missing = 0;

            if sfc_mode_group ~=5
                
                if ~do_partial_coherence
                    curr_trials = ctgsetli_good{j}(:,i);
                    Ntrials_min = ctgsetli_maxes{j}(i);
                
                    % Pull out chosen lfp traces
                    switch ue_pairs
                        case{0,2}
                            [currlfp1, lfpia_stage, adj_missing] = get_unitLFP_ts(lfp_sample,curr_unit,curr_stage,md,lfp_sample_indices,do_adjacent);
                            [currspike1] = get_spike_ts(curr_unit,curr_stage,md);
                            currlfp1_ctg0 = currlfp1(:,curr_trials);
                            currspike1_ctg0 = currspike1(:,curr_trials);
                            if plot_debug
                                plot_abs_vs_relative_lfp_and_spikes(filename_orig,lfp_sample_indices,curr_unit,md,currlfp1,currspike1,lfpia_stage);
                            end
                            clear currlfp1 currspike1
                        case 3; currspike1_ctg0 = currspike_all0(:,curr_trials,mypairs(j,1)); currlfp1_ctg0 = currlfp_all0(:,curr_trials,mypairs(j,2));
                        case 4; currlfp1_ctg0 = currlfp_all0(:,curr_trials,mypairs(j,1)); currlfp2_ctg0 = currlfp_all0(:,curr_trials,mypairs(j,2));
                        case 5; currspike1_ctg0 = currspike_all0(:,curr_trials,mypairs(j,1)); currspike2_ctg0 = currspike_all0(:,curr_trials,mypairs(j,2));
                        case 6; currlfp1_ctg0 = currlfp_all0(:,curr_trials,mypairs(j,1));
                            
                        case 7
                            [currspike1] = get_spike_ts(curr_unit,curr_stage,md); currspike1_ctg0 = currspike1(:,curr_trials);
                            clear currspike1
                    end
                    
                    if permutation_mode == 2
                        % Bootstrap a new distribution of trials             
                        Ntrials_curr = sum(curr_trials);
                        perm = floor(unifrnd(1,Ntrials_curr+1-0.0001,1,Ntrials_curr)); % Uniform sampling from 1 to Ntrials

                        % Test plotting...
                        % figure; plot(1:Ntrials_curr); hold on; plot(perm,'g'); hold on; plot(sort(perm),'r.'); legend('Original indices (identity permutation)','permuted','permuted sorted');

                        if exist('currspike1_ctg0','var'); currspike1_ctg0 = currspike1_ctg0(:,perm); end
                        if exist('currspike2_ctg0','var'); currspike2_ctg0 = currspike2_ctg0(:,perm); end
                        if exist('currlfp1_ctg0','var'); currlfp1_ctg0 = currlfp1_ctg0(:,perm); end
                        if exist('currlfp2_ctg0','var'); currlfp2_ctg0 = currlfp2_ctg0(:,perm); end
                    end
                end
                
            else
                switch mode_subgroups(1)    % Check if mode is to do unit analysis trace plotting
                    case 0
                        spikerates = get_numspikes_fromfile(curr_unit,curr_stage,md);
%                         % Normalizing
%                         numspikes = (numspikes - min(numspikes)) / (max(numspikes) - min(numspikes));
                        spikerates = spikerates(curr_trials);
                    case 1
                        [currspike1] = get_spike_ts(curr_unit,curr_stage,md);
                        currspike1 = currspike1(:,curr_trials);
                        spikerates = mean(currspike1,2); spikerates = spikerates / dt;
                        
                        apply_filter = 1;   % Apply sgolay filter to spiking data
                        if apply_filter
                            spikerates = sgolayfilt(spikerates,3,round(0.151/dt));   % Apply 3rd-order filter, 151 ms 
                        end
                end
            end
            
            % Basline subtraction
            if ~do_partial_coherence
                % If not doing partial, only baseline subtract the electrodes we're
                % considering
                switch ue_pairs
                    case {0,2,3};
                        currlfp1_ctg = do_baseline_subtract(currlfp1_ctg0,baseline_subtract);
                        currspike1_ctg = do_baseline_subtract(currspike1_ctg0,baseline_subtract);
                    case {4};
                        currlfp1_ctg = do_baseline_subtract(currlfp1_ctg0,baseline_subtract);
                        currlfp2_ctg = do_baseline_subtract(currlfp2_ctg0,baseline_subtract);
                    case {5};
                        currspike1_ctg = do_baseline_subtract(currspike1_ctg0,baseline_subtract);
                        currspike2_ctg = do_baseline_subtract(currspike2_ctg0,baseline_subtract);
                    case 6
                        currlfp1_ctg = do_baseline_subtract(currlfp1_ctg0,baseline_subtract);
                    case 7
                        currspike1_ctg = do_baseline_subtract(currspike1_ctg0,baseline_subtract);
                end
            end
            
            switch sfc_mode_group
                case 1
                    %[C ftemp] = sfc_fries2D(currspike(:,1),currlfp(:,1),dt,wind_s);
                    wind_s = round(0.2/dt);     % Window size for Fries method
                    
                    %tic
                    [Cave f lfp_wind] = sfc_fries2D(currspike1_ctg,currlfp1_ctg,dt,wind_s);
                    lfp_wind_mu = mean(lfp_wind,2);
                    lfp_wind_std = std(lfp_wind,[],2);
                    lfp_wind_ste = lfp_wind_std/sqrt(size(lfp_wind,2));
                    %toc
                    
                    sfc{i}.Cave = [sfc{i}.Cave Cave(:)];
                    sfc{i}.lfp_wind = [sfc{i}.lfp_wind lfp_wind_mu(:)];
                    sfc{1}.f = f;
                    
                    sfc{i}.S12ave = [];
                    sfc{i}.S1ave = [];
                    sfc{i}.S2ave = [];
                    sfc{i}.phiave = [];
                    
                    
                case 2
                    %tic
                    params.Fs = 1/dt;
                    params.trialave = 1;
                    
                    
                    
                    
                    if plot_debug1_thinning
                       [Cave2,phi,S12ave,S1ave,S2ave,f2,zerosp,confC,phistd,Cerr]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@coherencycpb,currlfp1_ctg,currspike1_ctg,params);
                       rate = mean(mean(currspike1_ctg0(:)))/dt;
                       [Cave3,phi,S12ave,S1ave,S2ave,f3,zerosp,confC,phistd,Cerr]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@Adjcoherencycpb,currlfp1_ctg,currspike1_ctg,params,[],[],StarRate,rate);
                     end
                    
                    %sfc{i}.spikes_per_trial{j} = sum(currspike_ctg,1)';
                    
                    %currlfp = zscore(currlfp); % Note- I found that taking the zscore doesn't make any difference. I guess this is because things are normalized already by the SFC operation.
                    
                    % Remove traces below a set number of spikes
%                     [currlfp_ctg,currspike_ctg] = trim_lowspike_ctgs(currlfp_ctg,currspike_ctg,thresh);
                    if thinning_mode == 0
                        [Cave,phi,S12ave,S1ave,S2ave,f,zerosp,confC,phistd,Cerr]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@coherencycpb,currlfp1_ctg,currspike1_ctg,params);
                    elseif thinning_mode == 1
                        rate = mean(mean(currspike1_ctg0(:)))/dt;
                        [Cave,phi,S12ave,S1ave,S2ave,f,zerosp,confC,phistd,Cerr]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@Adjcoherencycpb,currlfp1_ctg,currspike1_ctg,params,[],[],StarRate,rate);
                    elseif thinning_mode == 2
                        rate = mean(mean(currspike1_ctg0(:)))/dt;
                        prob_thinning = StarRate / rate;
                        if prob_thinning < 1;
                            [currlfp1_ctg, currspike1_ctg] = thin_spks_bootstrap(currlfp1_ctg, currspike1_ctg0, prob_thinning);     % Do thinning
                            currspike1_ctg = do_baseline_subtract(currspike1_ctg,baseline_subtract);                                % Recalculate baseline subtract, if necessary
                        end
                        [Cave,phi,S12ave,S1ave,S2ave,f,zerosp,confC,phistd,Cerr]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@coherencycpb,currlfp1_ctg,currspike1_ctg,params);
                    end
                    %toc
%                     figure('Position',[  1          26        1280         769]);
%                     range = 1
%                     subplot(221); title([filename ' unit # ' num2str(j) ' ' unit_names2])
%                     hold on; plot(f,Cave(:,range));legend('Cave'); xlabel('f (Hz)'); ylabel('SFC');
%                     subplot(222); hold on; plot(f,S1ave(:,range));legend('S1ave'); xlabel('f (Hz)'); ylabel('Power');
%                     subplot(223); hold on; plot(f,S2ave(:,range));legend('S2ave'); xlabel('f (Hz)'); ylabel('Power');
%                     subplot(224); hold on; plot(sum(currspike_ctg)); legend('Tot spikes in each trial'); xlabel('Trial #'); ylabel('# Spikes');
% 
% 
%                     figure; plot(f,Cave);
%                     
%                     figure
%                     for jj = 1:size(currlfp2_ctg,3)
%                         subplotsq(size(currlfp2_ctg,3),jj);
%                         title([num2str(jj) '-' convert_unit_underscores(md.lfp_names(jj))])
%                         [Cave3,phi2,S12ave2,S1ave2,S2ave2,f3]=coherencycpb(currlfp2_ctg(:,:,jj),currspike_ctg,params);
%                         hold on; plot(f3,Cave3); ylim([0 0.25]);
%                     end
%                     keyboard


                    if plot_debug1_thinning
                        clf; plot(f,Cave);
                        hold on; plot(f2,Cave2,'r');
                        hold on; plot(f3,Cave3,'g');
                        legend('Thinned','Default','Mikio');
                        pause(1);
                    end
                    
                    sfc{i}.Cave = [sfc{i}.Cave Cave(:)];
                    sfc{i}.S12ave = [sfc{i}.S12ave S12ave(:)];
                    sfc{i}.S1ave = [sfc{i}.S1ave S1ave(:)];
                    sfc{i}.S2ave = [sfc{i}.S2ave S2ave(:)];
                    sfc{i}.phiave = [sfc{i}.phiave phi(:)];
                    sfc{1}.f = f;
                    sfc{i}.confC = [sfc{i}.confC confC(:)];
                    sfc{i}.phistd = [sfc{i}.phistd phistd(:)];
                    
                    sfc{i}.Cerr1 = [sfc{i}.Cerr1 Cerr(1,:)'];
                    sfc{i}.Cerr2 = [sfc{i}.Cerr2 Cerr(2,:)'];
                    
                    
                    
%                     params.trialave = 0;
%                     [C,phi,S12,S1,S2,f]=coherencycpb(currlfp,currspike,params);
%                     C = C(:,sum(isnan(C)) == 0);
%                     Cave_cust = abs(mean(S12,2)./sqrt(mean(S1,2).*mean(S2,2)));
% 
%                     Aave = Cave; A=C; figure; plot(f,mean(Aave,2)); hold on; plot(f,mean(A,2),'r')
%                     xlabel('f (Hz)'); ylabel('Coherence');legend('Pre-average','Post-average');
%                     title(['F# '  filename 'mode ' num2str(sfc_mode_group) ' unit# ' num2str(curr_unit)])

                case 3
                    %tic
                    params.Fs = 1/dt;
                    params.trialave = 1;
                    
                    %sfc{i}.spikes_per_trial{j} = sum(currspike_ctg,1)';

                    % Remove traces below a set number of spikes
%                     [currlfp_ctg,currspike_ctg] = trim_lowspike_ctgs(currlfp_ctg,currspike_ctg,thresh);
                    if thinning_mode == 0
                        %[Ssp, phi, S12ave, S1ave, S2ave, F0, T, confC] = coherence_spect_wrapper(currlfp_ctg,currspike_ctg,'fs',1/dt,'params',params,'Nwind',Nwind);
                        
                        coherencycpb_Cerr = @(varargin) coherency_trials_downsample(Ntrials_min,baseline_subtract,@coherencycpb,varargin{:});
                        fname = @(X,Y) coherencycpb_Cerr(X,Y,params);
                        [Ssp, phi, S12ave, S1ave, S2ave, F0, zerosp, confC, phistd, Cerr, T] = func2spectrogram(currlfp1_ctg,currspike1_ctg,'fs',1/dt,'Nwind',Nwind,'fract_overlap',fract_overlap,'fname',fname);
                    elseif thinning_mode == 1
                        %[Ssp, phi, S12ave, S1ave, S2ave, F0, T, confC] = coherence_spect_wrapper(currlfp_ctg,currspike_ctg,'fs',1/dt,'params',params,'Nwind',Nwind,'mode',2.1);
                        
                        
                        coherencycpb_Cerr = @(varargin) coherency_trials_downsample(Ntrials_min,baseline_subtract,@Adjcoherencycpb,varargin{:});
                        fname = @(X,Y,Z) coherencycpb_Cerr(X,Y,params,[],[],StarRate,mean(Z(:))/dt);    % Z entry contains unnormalized spike rates.
                        [Ssp, phi, S12ave, S1ave, S2ave, F0, zerosp, confC, phistd, Cerr, T] = func2spectrogram(currlfp1_ctg,currspike1_ctg,currspike1_ctg0,'fs',1/dt,'Nwind',Nwind,'fract_overlap',fract_overlap,'fname',fname);
                        
                    elseif thinning_mode == 2
                        StarRate = get_rateStar;
                        rate = mean(mean(currspike1_ctg0(:)))/dt;
                        prob_thinning = StarRate / rate;
                        if prob_thinning < 1; [currlfp1_ctg, currspike1_ctg] = thin_spks_bootstrap(currlfp1_ctg, currspike1_ctg, prob_thinning); end
                        
                        coherencycpb_Cerr = @(varargin) coherency_trials_downsample(Ntrials_min,baseline_subtract,@coherencycpb,varargin{:});
                        fname = @(X,Y) coherencycpb_Cerr(X,Y,params);
                        [Ssp, phi, S12ave, S1ave, S2ave, F0, zerosp, confC, phistd, Cerr, T] = func2spectrogram(currlfp1_ctg,currspike1_ctg,'fs',1/dt,'Nwind',Nwind,'fract_overlap',fract_overlap,'fname',fname);
                        
                        %[Ssp, phi, S12ave, S1ave, S2ave, F0, T, confC] = coherence_spect_wrapper(currlfp_ctg,currspike_ctg,'fs',1/dt,'params',params,'Nwind',Nwind);
                        
                    end
                    F0 = mode(F0,2);
                    F0 = F0(:)'; T=T(:)';

                    
                    % Trim Ssp to only go up to 100 Hz
                    index = find(F0 <= 120);
                    F = F0(index);
                    Ssp = Ssp(index,:);
                    
                    % Trim Ssp to only go up to 250 Hz
                    index = find(F0 <= 250);
                    S1ave = S1ave(index,:);
                    S2ave = S2ave(index,:);
                    F2 = F0(index);

                    % Pack variables
                    Ssp = permute(Ssp,[1,3,2]);
                    S1ave = permute(S1ave,[1,3,2]);
                    S2ave = permute(S2ave,[1,3,2]);
                    confC = permute(confC,[1,3,2]);
                    
                    sfc{i}.Ssp = cat(2,sfc{i}.Ssp,Ssp);
                    %sfc{i}.S1AVE = cat(2,sfc{i}.S1ave,S1ave);    % Don't save this stuff to save space
                    sfc{i}.S2AVE = cat(2,sfc{i}.S2AVE,S2ave); 
                    sfc{i}.ConfC = cat(2,sfc{i}.ConfC,confC);

                    sfc{1}.F = F;
                    sfc{1}.F2 = F2;
                    sfc{1}.T = T;
                    
                case 4

                    params.Fs = 1/dt;
                    params.trialave = 0;
                    
                    
                    %currlfp = zscore(currlfp); % Note- I found that taking the zscore doesn't make any difference. I guess this is because things are normalized already by the SFC operation.
                    
                    % Remove traces below a set number of spikes
%                     [currlfp_ctg,currspike_ctg] = trim_lowspike_ctgs(currlfp_ctg,currspike_ctg,thresh);
                    [Cave,phi,S12ave,S1ave,S2ave,f]=coherencycpb(currlfp1_ctg,currspike1_ctg,params);
                    
                    %freqband = [20 30];
                    %index = (f >= freqband(1) & f < freqband(2));
                    %S12ave = convert_average(S12ave,index);
                    %S1ave = convert_average(S1ave,index);
                    %S2ave = convert_average(S2ave,index);
                    
                    index = find(f >= 25,1,'first');
                    S12ave = S12ave(index,:)';
                    S1ave = S1ave(index,:)';
                    S2ave = S2ave(index,:)';

                    sfc{i}.S12ave = [sfc{i}.S12ave S12ave(:)];
                    sfc{i}.S1ave = [sfc{i}.S1ave S1ave(:)];
                    sfc{i}.S2ave = [sfc{i}.S2ave S2ave(:)];
                    sfc{1}.f = f;
                    
                    
                    
%                     params.trialave = 0;
%                     [C,phi,S12,S1,S2,f]=coherencycpb(currlfp,currspike,params);
%                     C = C(:,sum(isnan(C)) == 0);
%                     Cave_cust = abs(mean(S12,2)./sqrt(mean(S1,2).*mean(S2,2)));
% 
%                     Aave = Cave; A=C; figure; plot(f,mean(Aave,2)); hold on; plot(f,mean(A,2),'r')
%                     xlabel('f (Hz)'); ylabel('Coherence');legend('Pre-average','Post-average');
%                     title(['F# '  filename 'mode ' num2str(sfc_mode_group) ' unit# ' num2str(curr_unit)])


                case 5
                    sfc{i}.spikerates = [sfc{i}.spikerates spikerates(:)];
                    
                case 6
                    [~, fbasename] = fileparts(filename_orig);
                    fsavepath = fullfile(getpath('path_buffer_curr'),['mode_' num2str(sfc_mode,12) fname_suffix]);
                    fsavename = fullfile(['f' num2str(filenum) 'u' num2str(j) '_' fbasename]);
                    
                    plot_sfc_on = 1;
                    plot_unit_summary(currspike1,curr_stage,ctgsetli,filename_orig,md.unit_names{j},curr_unit,plot_sfc_on,fsavepath,fsavename);
                    %plot_raster_with_sfc(currspike,curr_stage,ctgsetli,filename,md.unit_names{j});
                    
                case 8      % LFP Spectrogram
                    params.Fs = 1/dt;
                    params.trialave = 1;
                    [Ssp, F, T] = spect_wrapper(currlfp1_ctg,'fs',1/dt,'Nwind',Nwind,'params',params);
                    
                    % Trim Cave to only go up to 250 Hz
                    index = find(F <= 250);
                    F = F(index);
                    Ssp = Ssp(index,:);

                    % Pack variables
                    Ssp = permute(Ssp,[1,3,2]);
                    sfc{i}.Ssp = cat(2,sfc{i}.Ssp,Ssp);

                    sfc{1}.F = F;
                    sfc{1}.T = T;
                    
                case 9      % Evoked potential
                    E = mean(currlfp1_ctg,2);

                    % Pack variables
                    sfc{i}.E = cat(2,sfc{i}.E,E(:));
                    
                    
                case 12
                    
                    % =================================================================
                    % Parameters for Kyle modes
                    % =================================================================
                    bandwidth   = 10; % Hz.
                    start_f     = 10.0;     % Hz
                    stop_f      = 100.0;    % Hz
                    n_trim      = 100;
                    f_out       = 'out/glm_hs_simple.out'; 
                    debug_level = 0;  % To get phase output set >= 1.

                    v2  = m_form_v( currspike1_ctg,currlfp1_ctg);
                    
                    
                    % ========================================
                    % GBC -- piecewise linear
                    % ========================================
                    h_func        = @(x) ( max( x, 0 ));
                    h_func_inv    = @(x) ( max( x, 0 ));
                    h_func_deriv  = @(x) ( ones(size(x)) );

                    
                    gbc_out_pl  = gbc( v2, dt,                              ...
                                     h_func, h_func_deriv, h_func_inv,    ...
                                     bandwidth, start_f, stop_f, n_trim,  ...
                                     f_out, debug_level );
                    
                    sout = gbc_out_pl;
                    
                    % Old fields (same as Chronux)
                    sfc{i}.Cave = [sfc{i}.Cave sout.rho(:)];
                    sfc{i}.S12ave = [sfc{i}.S12ave []];
                    sfc{i}.S1ave = [sfc{i}.S1ave []];
                    sfc{i}.S2ave = [sfc{i}.S2ave []];
                    sfc{i}.phiave = [sfc{i}.phiave sout.phi_p_in_deg(:)/180*pi];        % Phi should be in radians (-pi to pi)
                    sfc{1}.f = sout.f;
                    % New fields
                    sfc{i}.alpha = [];
                    
                    
                    
                case 15
                    % =================================================================
                    % Parameters for Kyle modes
                    % =================================================================
                    bandwidth   = 10; % Hz.
                    start_f     = 20.0;     % Hz
                    stop_f      = 100.0;    % Hz
                    n_trim      = 100;
                    f_out       = 'out/glm_hs_simple.out'; 
                    debug_level = 0;  % To get phase output set >= 1.

                    v2  = m_form_v( currspike1_ctg,currlfp1_ctg);


                    log_out     = glm_sfa_log_dave( v2, dt, ...
                                           bandwidth, start_f, stop_f, n_trim, ...
                                           f_out, debug_level );

                    sout = log_out;

                    rho_log  = sqrt( log_out.coeffs(2,:).^2 + log_out.coeffs(3,:).^2 );
                    alpha_log = log_out.coeffs(1,:);
                    phipref = atan2( log_out.coeffs(3,:),log_out.coeffs(2,:));        % Phi should be in radians (-pi to pi)

                    
                    % Old fields (same as Chronux)
                    sfc{i}.Cave = [sfc{i}.Cave rho_log(:)];
                    sfc{i}.S12ave = [sfc{i}.S12ave []];
                    sfc{i}.S1ave = [sfc{i}.S1ave []];
                    sfc{i}.S2ave = [sfc{i}.S2ave []];
                    sfc{i}.phiave = [sfc{i}.phiave phipref(:)];
                    sfc{i}.p = [sfc{i}.p sout.p_zeromod(:)];
                    sfc{1}.f = sout.f;
                    % New fields
                    sfc{i}.alpha = [sfc{i}.alpha alpha_log(:)];
                    
                
                case 22         % FFC
                    params.Fs = 1/dt;
                    params.trialave = 1;
                    
                    if ~do_partial_coherence
                        [Cave,phi,S12ave,S1ave,S2ave,f,confC,phistd,Cerr]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@coherencyc,currlfp1_ctg,currlfp2_ctg,params);
                    else
                        %params.fname = @coherencyc;
                        %params.chosen_dims = [mypairs(j,1) mypairs(j,2)];
                        %[Cave,phi,S12ave,S1ave,S2ave,f,confC,phistd,Cerr]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@partcoh_stupidwrapper,currlfp_all_ctg,currlfp_all_ctg,params);
                        Cave = Cave_all(:,j); phi = phiall(:,j);
                        confC=[];phistd=[];Cerr=[NaN;NaN];
                        S1ave=[];S2ave=[];S12ave=[];
                        
                        ctgsetli_good{j}(:,i) = ctgsetli_good_all(:,i);
                        
                        
                    end
                    
                    sfc{i}.Cave = [sfc{i}.Cave Cave(:)];
                    sfc{i}.S12ave = [sfc{i}.S12ave S12ave(:)];
                    sfc{i}.S1ave = [sfc{i}.S1ave S1ave(:)];
                    sfc{i}.S2ave = [sfc{i}.S2ave S2ave(:)];
                    sfc{i}.phiave = [sfc{i}.phiave phi(:)];
                    sfc{1}.f = f;
                    sfc{i}.confC = [sfc{i}.confC confC(:)];
                    sfc{i}.phistd = [sfc{i}.phistd phistd(:)];
                    
                    sfc{i}.Cerr1 = [sfc{i}.Cerr1 Cerr(1,:)'];
                    sfc{i}.Cerr2 = [sfc{i}.Cerr2 Cerr(2,:)'];
                    
                    
                case 23
                    params.Fs = 1/dt;
                    params.trialave = 1;
                    
                    if ~do_partial_coherence
                        %[Cave,phi,S12ave,S1ave,S2ave,f,confC,phistd,Cerr]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@coherencyc,currlfp1_ctg,currlfp2_ctg,params);
                        
                        coherencyc_Cerr = @(varargin) coherency_trials_downsample(Ntrials_min,baseline_subtract,@coherencyc,varargin{:});
                        fname = @(X,Y) coherencyc_Cerr(X,Y,params);
                        [Ssp, phi, S12ave, S1ave, S2ave, F, confC, phistd, Cerr, T] = func2spectrogram(currlfp1_ctg,currlfp2_ctg,'fs',1/dt,'Nwind',Nwind,'fract_overlap',fract_overlap,'fname',fname);
                        
                    else
                        % Never implemented
                    end
                    
                    f = F(:,1);
                    
                    % Trim Ssp to only go up to 100 Hz
                    index = find(f <= 120);
                    f = f(index);
                    Ssp = Ssp(index,:);
                    phi = phi(index,:);
                    phistd = phistd(index,:);
                    
                    sfc{i}.Cave = cat(2,sfc{i}.Cave,Ssp(:));     % Spectrum data
                    sfc{i}.phiave = [sfc{i}.phiave phi(:)];
                    sfc{1}.f = f;
                    sfc{1}.f2 = T;
                    sfc{i}.confC = [sfc{i}.confC confC(1)];
                    sfc{i}.phistd = [sfc{i}.phistd phistd(:)];
                    
                case 32         % SSC

                    params.Fs = 1/dt;
                    params.trialave = 1;
                    
                    if plot_debug1_thinning
                         [Cave2,phi,S12ave,S1ave,S2ave,f2,zerosp,confC,phistd,Cerr]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@coherencypb,currspike1_ctg,currspike2_ctg,params); 

                         rate1 = mean(mean(currspike1_ctg0(:)))/dt; 
                         rate2 = mean(mean(currspike2_ctg0(:)))/dt; 
                         [Cave3,phi,S12ave,S1ave,S2ave,f3,zerosp,confC,phistd,Cerr]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@Adjcoherencypb,currspike1_ctg,currspike2_ctg,params,[],StarRate,rate1,StarRate,rate2);  % NEeds updating to Mikio's method
                     end
                    
                    % To be used if I want to look at correlation sometime
%                     spikerates1 = mean(currspike1_ctg0,2); spikerates1 = spikerates1 / dt;
%                     spikerates1 = sgolayfilt(spikerates1,3,round(0.151/dt));   % Apply 3rd-order filter, 151 ms 
%                     spikerates2 = mean(currspike2_ctg0,2); spikerates2 = spikerates2 / dt;
%                     spikerates2 = sgolayfilt(spikerates2,3,round(0.151/dt));   % Apply 3rd-order filter, 151 ms 
                                     
                    if thinning_mode == 0
                        [Cave,phi,S12ave,S1ave,S2ave,f,zerosp,confC,phistd,Cerr]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@coherencypb,currspike1_ctg,currspike2_ctg,params); 
                    elseif thinning_mode == 1
                        rate1 = mean(mean(currspike1_ctg0(:)))/dt; 
                        rate2 = mean(mean(currspike2_ctg0(:)))/dt; 
                        [Cave,phi,S12ave,S1ave,S2ave,f,zerosp,confC,phistd,Cerr]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@Adjcoherencypb,currspike1_ctg,currspike2_ctg,params,[],StarRate,rate1,StarRate,rate2);  % NEeds updating to Mikio's method
                    elseif thinning_mode == 2
                        
                        % % % % % % % % This is not working! % % % % % % %
                        rate1 = mean(mean(currspike1_ctg0(:)))/dt; prob_thinning1 = StarRate / rate1;
                        rate2 = mean(mean(currspike2_ctg0(:)))/dt; prob_thinning2 = StarRate / rate2;
                        
                        sz = size(currspike2_ctg0);
                        [Nbs1] = est_Nbootstraps(prob_thinning1, sz(2));
                        [Nbs2] = est_Nbootstraps(prob_thinning2, sz(2));
                        Nbootstraps = min(Nbs1,Nbs2);
                        clear Nbs1 Nbs2
                        
                        if Nbootstraps > 0
                            [~, currspike1_ctg] = thin_spks_bootstrap([], currspike1_ctg0, prob_thinning1,Nbootstraps);     % Do thinning
                            [~, currspike2_ctg] = thin_spks_bootstrap([], currspike2_ctg0, prob_thinning2,Nbootstraps);     
                            
                            currspike1_ctg = do_baseline_subtract(currspike1_ctg,baseline_subtract);           % Recalculate baseline subtract, if necessary
                            currspike2_ctg = do_baseline_subtract(currspike2_ctg,baseline_subtract);           
                        end
                        
                        [Cave,phi,S12ave,S1ave,S2ave,f,zerosp,confC,phistd,Cerr]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@coherencypb,currspike1_ctg,currspike2_ctg,params);
                        
                    end
                    
                    if plot_debug1_thinning
                        clf; plot(f,Cave);
                        hold on; plot(f2,Cave2,'r');
                        hold on; plot(f3,Cave3,'g');
                        pause(1);
                    end
                    
                    sfc{i}.Cave = [sfc{i}.Cave Cave(:)];
                    sfc{i}.S12ave = [sfc{i}.S12ave S12ave(:)];
                    sfc{i}.S1ave = [sfc{i}.S1ave S1ave(:)];
                    sfc{i}.S2ave = [sfc{i}.S2ave S2ave(:)];
                    sfc{i}.phiave = [sfc{i}.phiave phi(:)];
                    sfc{1}.f = f;
                    sfc{i}.confC = [sfc{i}.confC confC(:)];
                    sfc{i}.phistd = [sfc{i}.phistd phistd(:)];
                    
                    sfc{i}.Cerr1 = [sfc{i}.Cerr1 Cerr(1,:)'];
                    sfc{i}.Cerr2 = [sfc{i}.Cerr2 Cerr(2,:)'];
                    
                case 40                 % This calculates the cross frequency coupling
                    params.Fs = 1/dt;
                    currlfp1_ctg;
                    params2 = params;
                    params2.trialave = 0;
                    [S,f,Serr] = mtspectrumc(currlfp1_ctg,params2);
                    
                    % Do CFC
                    ind = (f >= 2) & (f <= 65);
                    S=S(ind,:)';
                    f=f(ind);
                    [R,P,RLO,RUP] = corrcoef(S);
                    
                    
%                     figure; imagesc(R)
%                     figure; imagesc(f,f,(P < 0.05 & R < 0).*R); set(gca,'YDir','normal')
                    
                    R = permute(R,[1,3,2]);
                    P = permute(P,[1,3,2]);
                    RLO = permute(RLO,[1,3,2]);
                    RUP = permute(RUP,[1,3,2]);
                    sfc{i}.Cave = cat(2,sfc{i}.Cave,R);      % Use Cave label for our primary output
                    sfc{i}.Ptr = cat(2,sfc{i}.Ptr,P);
                    sfc{i}.Cerr1 = cat(2,sfc{i}.Cerr1,RLO);
                    sfc{i}.Cerr2 = cat(2,sfc{i}.Cerr2,RUP);
                    sfc{1}.f = f;
                    
                    

                case 41                 % Just calculate power spectrum
                    params.Fs = 1/dt;
                    currlfp1_ctg;
                    params2 = params;
                    params2.trialave = 1;
                    [S,f,Serr] = mtspectrumc(currlfp1_ctg,params2);

                    
%                     figure; imagesc(R)
%                     figure; imagesc(f,f,(P < 0.05 & R < 0).*R); set(gca,'YDir','normal')
                    
                    sfc{i}.Cave = cat(2,sfc{i}.Cave,S);     % Spectrum data
                    sfc{i}.Cerr1 = cat(2,sfc{i}.Cerr1,Serr(1,:)');
                    sfc{i}.Cerr2 = cat(2,sfc{i}.Cerr2,Serr(2,:)');
                    sfc{1}.f = f;
                    
                case {42,43,44}
                    % Necessary inputs values
                    sig_pac = currlfp1_ctg(:,:);
                    Fs = 1/dt;
                    switch sfc_mode_group
                        case 42
                            measure = 'esc';
                        case 43
                            measure = 'mi';
                        case 44
                            measure = 'cfc';
                    end
                            
                    % Other input parameters (mostly defaults)
                    sig_mod = sig_pac;
                   ph_freq_vec = 2:3:29;    % Full range
                   amp_freq_vec = 2:5:101; 
%                     ph_freq_vec = 2:1:29;    % Expanded range
%                     amp_freq_vec = 10:3:160; 
%                     ph_freq_vec = 10:3:20;   % Sub range
%                     amp_freq_vec = 15:10:55;
%                     ph_freq_vec = [18 20]; % Single value
%                     amp_freq_vec = [40 45];
                    plt = 'n';
                    waitbar = 0;
                    width = 7;
                    nfft = ceil(Fs/(diff(ph_freq_vec(1:2))));
                    if tapers_mode == 2
                        num_shf = 150;
                    else
                        num_shf = 0;
                    end
                    alpha = 0.05;
                    dataname = '';
                    sig_pac_name = '';
                    sig_mod_name = '';

                    % Calculate PAC
                    [pacmat, freqvec_ph, freqvec_amp, p] = find_pac_shf (sig_pac, Fs, measure, ...
                    sig_mod, ph_freq_vec, amp_freq_vec, plt, waitbar, width, nfft, num_shf, alpha,...
                    dataname, sig_pac_name, sig_mod_name);
                
                    
                    % Replace upper triangle with NaNs to save space. Only
                    % works with square matrices
%                     sub = upper_triangle_pairs(size(pacmat,1));
%                     ind = sub2ind(size(pacmat),sub(:,1),sub(:,2));
%                     pacmat_orig=pacmat; p_orig = p;
%                     pacmat(ind) = NaN;
%                     p(ind) = NaN;
                    
%                     figure; imagesc(freqvec_ph,freqvec_amp,pacmat); caxis([0 500]); colorbar

                    sfc{i}.Cave = cat(2,sfc{i}.Cave,pacmat(:));     % Spectrum data
                    sfc{i}.p = cat(2,sfc{i}.p,p(:));     % Spectrum data
                    sfc{1}.f = ph_freq_vec;
                    sfc{1}.f2 = amp_freq_vec;
                    
                case 52
                    switch tapers_mode
                        case 0          % Average entire spike train
                            
                            % % More complicated way of doing things % %
                            % spikerates2 = currspike1_ctg / dt; % Converts spike event (0 or 1) into spike rate per time bin.
                            % spikerates2 = mean(spikerates2,2); % Mean spike rate across trials
                            % spikerates2 = mean(spikerates2(:)); % Mean spike rate across all time.
                            
                            % % % This should be equal to the above % % %
                            spikerates2 = mean(sum(currspike1_ctg,1) / (size(currspike1_ctg,1)*dt)); % Should work out to the same as this.
                            sfc{i}.Cave = [sfc{i}.Cave repmat(spikerates2,5,1)]; 
                            sfc{1}.f = 1:5;
                        case {1,2,3,4}
                            Nwind = os.Nwind;
                            
                            spikerates2 = mean(currspike1_ctg,2) / dt;
                            spikerates2 = sgolayfilt(spikerates2,3,Nwind);   % Apply 3rd-order filter, Nwind samples long!
                            sfc{i}.Cave = [sfc{i}.Cave spikerates2];
                            sfc{1}.f = [1:length(spikerates2)]*dt;
                    end
                        
                    
            end
                    
                    
            
            if sfc_mode_group ~=5
                sfc{i}.Ntraces = [sfc{i}.Ntraces sum(curr_trials)];
                sfc{i}.mypairs = [sfc{i}.mypairs mypairs(j,:)'];
                %sfc{i}.sum_ctgsetli2 = [sfc{i}.sum_ctgsetli2 sum(ctgsetli)'];   % Actually this isn't necessary; we have Ntraces to save this info
                
                sfc{i}.sum_ctgsetli_good = [sfc{i}.sum_ctgsetli_good sum(ctgsetli_good{j}(:,i))];
                %sfc{i}.sum_ctgsetli_maxes = [sfc{i}.sum_ctgsetli_maxes ctgsetli_maxes{j}(i)];
                sfc{i}.sum_ctgsetli_maxes = [sfc{i}.sum_ctgsetli_maxes Ntrials_min];
                
                if ue_pairs ~= 4 && ue_pairs ~= 6           % Only save spike info if we're not in LFP only mode
                    sfc{i}.Ntraces_below_thresh = [sfc{i}.Ntraces_below_thresh sum(sum(currspike1_ctg0,1) < thresh)];
                    sfc{i}.spikerate_mu = [sfc{i}.spikerate_mu sum(sum(currspike1_ctg0,1)) / (size(currspike1_ctg0,1)*get_dt*sfc{i}.Ntraces(end))];    % total # spikes / total time
                end
                
                if ue_pairs == 5
                    sfc{i}.Ntraces_below_thresh2 = [sfc{i}.Ntraces_below_thresh2 sum(sum(currspike2_ctg0,1) < thresh)];
                    sfc{i}.spikerate_mu2 = [sfc{i}.spikerate_mu2 sum(sum(currspike2_ctg0,1)) / (size(currspike2_ctg0,1)*get_dt*sfc{i}.Ntraces(end))];    % total # spikes / total time
                end
                
                sfc{i}.adj_was_missing = [sfc{i}.adj_was_missing, adj_missing];
                sfc{i}.elapsedTime = [sfc{i}.elapsedTime, toc(ticID)];
            end
            
        end
    end
    
    if sfc_mode_group ~= 6
        if Ncoherences >= min(j0)        % Make sure there is at least 1 unit in the range of j specitied
            % Save output, unless we're in raster plotting mode.
            sfc{1}.singletonfieldnames = {'singletonfieldnames','f','f2' ...
                'F','F2','T', ...                                                   % For 3D plots (coherograms)
                'unit_names','filename','tapers','sum_ctgsetli',...
                'include_switch_trials','N_ctgs_base','N_ctgs_extras','Nwind','stagesir', ...
                'md_lfp_names','md_unit_names',...
                'spikethresh','StarRate','params','i0','j0'};
            sfc{1}.unit_names = unit_names;
            sfc{1}.filename = filename_orig;
            sfc{1}.tapers = get_tapers;
            sfc{1}.sum_ctgsetli = sum(ctgsetli0);
            sfc{1}.include_switch_trials = ctgsetli_mode; % Metadata
            sfc{1}.N_ctgs_base = N_ctgs_base;
            sfc{1}.N_ctgs_extras = N_ctgs_extras;
            sfc{1}.Nwind = Nwind;
            sfc{1}.stagesir = get_stagesir(curr_stage);
            
            % Original names
            sfc{1}.md_lfp_names = md.lfp_names;
            sfc{1}.md_unit_names = md.unit_names;

            % More stuff
            sfc{1}.spikethresh = thresh;
            sfc{1}.StarRate = StarRate;
            sfc{1}.params = params;
            sfc{1}.i0 = i0;
            sfc{1}.j0 = j0;

        end
    end

end





function plot_abs_vs_relative_lfp_and_spikes(filename,lfp_sample_indices,curr_unit,md,currlfp,currspike,lfpia_stage)
    %% function plot_abs_vs_relative_lfp_and_spikes(filename,lfp_sample_indices,curr_unit,md,currlfp,currspike,lfpia_stage)

    % % Plot everything
    % Plot lfp
    curr_electrode = md.unit_LFP_pairs(2,curr_unit);
    
    dt = get_dt;
    s=load (fullfile(getpath('path_lfp'),filename));     % Load full file
    lfpmat2 = double(s.lfpmat_int16); clear s;

    shorten_factor=50;

    lfpmat3 = lfpmat2(1:round(end/shorten_factor),curr_electrode);
    t2 = (1:length(lfpmat3))*dt;
    %t2 = (0:length(lfpmat3)-1)*dt;

    shorten_factor=100;
    lfpia = recall_lfpia(lfp_sample_indices);
    lfpia = lfpia(:,1:round(end/shorten_factor));

    figure; h{1}=plot(t2,lfpmat3);
    hold on; h{2}=plot(t2(lfpia),lfpmat3(lfpia),'r');



    % Overlay the units
    unitstemp = md.unit_ind{curr_unit};
    unitstemp = unitstemp(1:round(end/shorten_factor));
%     hold on; h{3}=plot(t2(unitstemp),500*ones(1,length(unitstemp)),'ro');    % Overlay the units at a fixed verticla position
    hold on; plot(t2(unitstemp),lfpmat3(unitstemp),'ko');               % Overlaythe units ontop of the LFP. Harder to read...


    % % Plot just sample on stuff
    % LFP
    currlfp2 = currlfp(:,1:round(end/shorten_factor));
    currspike2 = currspike(:,1:round(end/shorten_factor));
    lfpia_stage2 = lfpia_stage(:,1:round(end/shorten_factor));
    hold on; h{4}=plot(t2(lfpia_stage2),currlfp2,'g');
%     hold on; h{5}=plot(t2(lfpia_stage2),500*currspike2,'g.');
    
    
    % % Plot sample on times
    
    hold on;
    sot = md.sample_on_times(1:round(end/shorten_factor));
    h{6}=plot(sot,zeros(1,length(sot)),'kx');
    soi = md.sample_on_ind(1:round(end/shorten_factor));
    h{7}=plot(t2(soi),zeros(1,length(sot)),'bx');
    
% 
%     params.Fs = 1/dt;
% 
%     params.trialave = 1;
%     [Cave,phi,S12ave,S1ave,S2ave,f]=coherencycpb(currlfp2,currspike2,params);
%     Aave = Cave;figure; plot(f,mean(Aave,2));
%     xlabel('f (Hz)'); ylabel('Coherence');legend('Pre-average','Post-average');
%     title(['F# '  filename ' unit# ' num2str(curr_unit)])
%     
%     
%     [C f] = sfc_fries([currspike(:) currlfp(:)],dt,round(0.5/dt));
    

    clear lfpmat2 t2 lfpir lfpia lfp_indices_abs unitstemp currlfp2
end





function [currlfp_ctg,currspike_ctg] = trim_lowspike_ctgs(currlfp_ctg,currspike_ctg,thresh);
    %% function [currlfp_ctg,currspike_ctg] = trim_lowspike_ctgs(currlfp_ctg,currspike_ctg,thresh);
    index = sum(currspike_ctg) >= thresh;
    currlfp_ctg = currlfp_ctg(:,index);
    currspike_ctg = currspike_ctg(:,index);
end




function curr_badtrials = get_curr_badtrials(bad_trials,filenum,curr_unit,md)
    %% function curr_badtrials = get_curr_badtrials(bad_trials,filenum,curr_unit,md)
    curr_lfp = md.unit_LFP_pairs(2,curr_unit);
    curr_badtrials = bad_trials{filenum}(:,curr_lfp);
end




function X = convert_average(X,index)
    %% function X = convert_average(X,index)
    % Converts X into an average in a specific frequency band
    X = X(index,:);
    X = mean(X);
    X = X(:);
end







function ctgsetli_new = permute_ctgsetli(ctgsetli, i0)
    %% function ctgsetli_new = permute_ctgsetli(ctgsetli, i0)
    test_on = 0;

    % These are the two groups prior to permutation
    ctgs1 = ctgsetli(:,i0(1)); N1 = sum(ctgs1);
    ctgs2 = ctgsetli(:,i0(2)); N2 = sum(ctgs2);
    
    % This is the pool, from which we will be drawing
    samples_pool = ctgs1 | ctgs2; Nall = sum(samples_pool);
    
    ind = find(samples_pool);
    ind = ind(randperm(length(ind)));
    chosen1 = ind(1:N1);
    chosen2 = ind(N1+1:end);
    
    ctgsetli_new = ctgsetli;
    ctgsetli_new(:,i0(1)) = false;
    ctgsetli_new(:,i0(2)) = false;
    
    ctgsetli_new(chosen1,i0(1)) = true;
    ctgsetli_new(chosen2,i0(2)) = true;
    
    if test_on
    
        figure; imagesc(ctgsetli(:,i0)); title('Original');
        figure; imagesc(ctgsetli_new(:,i0)); title('New, after permutation')
        figure; plot(any(ctgsetli(:,i0),2),'bo'); hold on; plot(any(ctgsetli_new(:,i0),2),'r.'); title('New vs old samples pool (the "or" across rows of prev 2 figures; should be same'); legend('Ctgset','Ctgset new');
        
        sum(ctgsetli)
        sum(ctgsetli_new)
        
    end

end





function currspike_all = get_spike_ts_all(Nunits,curr_stage,md)
    %% function currspike_all = get_spike_ts_all(Nunits,curr_stage,md)

    j = 1;
    [currspike] = get_spike_ts(j,curr_stage,md);
    sz = size(currspike);
    currspike_all = zeros([sz, Nunits]);
    currspike_all(:,:,j) = currspike;
    
    for j = 1:Nunits
        currspike_all(:,:,j) = get_spike_ts(j,curr_stage,md);
    end
    
end







function out6 = calc_ctgsetli_mins_all(ctgsetli_good_A)
    %% function out6 = calc_ctgsetli_mins_all(ctgsetli_good_A)

%     % Method 1
%     tic
%     out1 = cellfunu(@(x) min(vertcat(x{:}),[],1),ctgsetli_good_A);   % Min for each file
%     out2 = vertcat(out1{:});       % Min for each file stacked in matrix
%     out3 = min(out2,[],1);
%     toc
    
    % Method 2 (this one seems a bit faster)
%     tic
    for i = 1:length(ctgsetli_good_A)
        out4{i} = vertcat(ctgsetli_good_A{i}{:});
    end
    out5 = vertcat(out4{:});
    out6 = min(out5,[],1);
%     toc
    
%     sum(out3 ~= out6);

end

function ctgs_out = ctgsetli_randtrim(ctgs,maxes)
    %% function ctgs_out = ctgsetli_randtrim(ctgs,maxes)
    plot_test=0;
    Nctgs = size(ctgs,2);
    ctgs_out = false(size(ctgs));
    
    for i  = 1:Nctgs
        inds=find(ctgs(:,i));
        inds = inds(randperm(length(inds)));    % Shuffle inds
        inds = inds(1:maxes(i));
        ctgs_out(inds,i) = true;
    end
    
    if plot_test
        figure
        subplot(121);imagesc(ctgs); title('Original Ctgsetli'); xlabel('Ctgs'); ylabel('Trials'); subplot(122);bar(sum(ctgs)); xlabel('Ctgs'); ylabel('Sum(trials)')
        figure;subplot(121);imagesc(ctgs_out); title('Trimmed Ctgsetli'); xlabel('Ctgs'); ylabel('Trials');
        subplot(122);bar([sum(ctgs_out);maxes]'); ylabel('Sum(trials)')
    end
end



function [ctgsetli_good, ctgsetli_maxes] = calc_ctgsetli_goods(j0_precalc,Ncoherences,sfc_mode_group,exclude_clipping_trials,bad_trials,filenum,md,ctgsetli0,ue_pairs,permutation_mode,mypairs,ctgsall,coh_debias_mode,i0)
    %% function [ctgsetli_good, ctgsetli_maxes] = calc_ctgsetli_goods(j0_precalc,Ncoherences,sfc_mode_group,exclude_clipping_trials,bad_trials,filenum,md,ctgsetli0,ue_pairs,permutation_mode,mypairs,ctgsall,coh_debias_mode,i0)
    plot_debug3_currtrials = 0;
    for j = j0_precalc
        if j > Ncoherences; continue; end
        curr_unit = j;

        if sfc_mode_group ~=5

            % Identify clipping trials
            if exclude_clipping_trials
                switch ue_pairs
                    case 0
                        curr_lfp = md.unit_LFP_pairs(2,curr_unit);
                        curr_badtrials = bad_trials{filenum}(:,curr_lfp);
                    case 2
                        [curr_lfp]= get_adjacent_electrode(md,curr_unit);
                        curr_badtrials = bad_trials{filenum}(:,curr_lfp);
                    case 3
                        curr_lfp = mypairs(j,2);
                        curr_badtrials = bad_trials{filenum}(:,curr_lfp);
                    case 4
                        curr_lfp = mypairs(j,1);
                        curr_lfp2 = mypairs(j,2);
                        curr_badtrials = bad_trials{filenum}(:,curr_lfp) | bad_trials{filenum}(:,curr_lfp2);
                    case {5,7}
                        %Nsamples = size(bad_trials{filenum},1);
                        Nsamples = size(ctgsetli0,1);
                        curr_badtrials = false(Nsamples,1); % All trials are good!
                    case 6
                        curr_lfp = j;
                        curr_badtrials = bad_trials{filenum}(:,curr_lfp);
                end
            else
                curr_badtrials = false(Nsamples,1); % All trials are good!
            end

            if plot_debug3_currtrials
                temp = currlfp_all0(:,:,curr_lfp);
                temp2=  bad_trials{filenum}(:,curr_lfp);
                n = 1:length(temp(:));
                n=reshape(n,size(temp));
                figure; plot(n,temp,'b');
                hold on; plot(n(:,temp2),temp(:,temp2),'r:');
            end


            ctgsetli_good{j} = ctgsetli0 & repmat(~curr_badtrials(:),1,size(ctgsetli0,2));


            % Calculate 
            if ~isempty(ctgsall)
                if any(ctgsall.ctgsetli_good_A{filenum}{j} ~= sum(ctgsetli_good{j}))
                    warning('Ctgsetli calculation mismatch');
                end
            end
            ctgsetli_maxes{j} = calc_ctgsetli_maxes(ctgsetli_good{j},ctgsall,coh_debias_mode);
            %ctgsetli_good{j} = ctgsetli_randtrim(ctgsetli_good{j},ctgsetli_maxes{j});


            if permutation_mode == 1
                if length(i0) ~= 2; warning('When permutation_mode=1, length(i0) must be 2.'); return; end
                ctgsetli_good{j} = permute_ctgsetli(ctgsetli_good{j}, i0);
            end
        end
    end
end
