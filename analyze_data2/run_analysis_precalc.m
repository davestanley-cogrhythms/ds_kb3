
function run_analysis_precalc (sfc_mode, curr_stage, file_range,i0,j0,outname_embed_i0j0)
%% function run_analysis_precalc (sfc_mode, stage, file_range,i0,j0,outname_embed_i0j0)

format compact
format long g

run_setdefaultfig
if ~isdeployed
    %% load paths
    addpath(genpath('./funcs_supporting_local'));
    addpath(genpath('./funcs_run_analysis_precalc'));
end

% if ~exist('sfc_mode','var'); sfc_mode = 22.401411101; end
% if ~exist('sfc_mode','var'); sfc_mode = 52.700201001; end
if ~exist('sfc_mode','var'); sfc_mode = 45.6018101041; end
% if ~exist('sfc_mode','var'); sfc_mode = 45.601811101; end
% if ~exist('sfc_mode','var'); sfc_mode = 52.70020100; end
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
    % sfc_mode = 23 -> LFP using Chronux Spectrogram
    
    % 30 <= sfc_mode < 40 - SFC all elects and FFC comparisons
    % sfc_mode = 32 -> SSC (all elects) 
    % sfc_mode = 33 -> SSC spectrogram (not implemented!)
    
    % sfc_mode = 40 -> Amplitude-Amplitude coupling, trial by trial
    % sfc_mode = 41 -> PSD amplitude
    % sfc_mode = 42 -> Cross-frequency coupling ESC method
    % sfc_mode = 43 -> Cross-frequency coupling MI method
    % sfc_mode = 44 -> Cross-frequency coupling CFC method
    % sfc_mode = 45 -> PSD spectrogram (not implemented!)
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
    %     X=5 analyse boundary trials SchA/B;
    %     X=6 boundary trials SchA 50% vs SchB 50% (no 100% trials);
    %     X=7 As mode0, but breaks down into 60% and 100% ;
    %     X=8 ??? I forget;
    %     X=9 For RoyFig (all morph percentages);
    %     X=10 As 0, but merge rel/irrel;
    %     X=40 As 0, but pair rel/irrel (so rel and irrel are paired: ctg1, ctg1 irrel, ctg2, ctg2 irrel, etc)
    %     X=24 analyse boundary trials (balanced);
    %     X=34 analyse boundary trials (balanced + merged rel/irrel);
    %     X=25 analyse boundary trials Sch A/B (balanced);
    %     X=35 analyse boundary trials Sch A/B (balanced + merged rel/irrel);
    %     
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
    
if ~exist('curr_stage','var'); curr_stage = 4; end
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
if ~exist('i0','var'); i0 = [1:8]; end       % There are a max of 11 categories Roy data; There are fewer for Cromer
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
bad_trials = prep_badtrials(curr_stage);

if  ~exist('file_range','var')
    file_range = 1:Nfiles;
%    file_range = [4];
%    file_range = [44];
%    file_range = [70];
%     file_range = [1];
%     file_range = [1:Nfiles];
end

for filenum = file_range
    
    %profile on
    
    % File name
    fname = file_list{filenum};
    fprintf ('Processing file %d of %d, name %s \n',filenum,file_range(end),fname);
    
    % Load options structure
    os = get_options(sfc_mode, curr_stage,i0,j0,outname_embed_i0j0,file_list,filenum);
    
    % Check if output name exists
    outpath = os.outpath;
    outname = os.outname;
    mkdir_dave(outpath);
    if exist(os.outname,'file')
        fprintf('File %s exists...skipping\n',outname);
        continue
    end
    
    % Load metadata
    md = load (fullfile(getpath('path_metadata'),fname));     % Load metadata
    
    % Load data (& tile if doing spectrograms)
    ticID = tic;
    if ~os.is_spectrogram
        [spike_all0,lfp_all0] = preload_and_save(os,md);
        do_spikes = 1; [J0_spike,Jparams_spike] = precalc_fft_and_save(spike_all0,os,do_spikes);
        do_spikes = 0; [J0_lfp,Jparams_lfp] = precalc_fft_and_save(lfp_all0,os,do_spikes);
        if os.needed_spikes; Jparams = Jparams_spike; elseif os.needed_lfp; Jparams = Jparams_lfp; end        % Should at least pick one of these.
        
        if os.needed_spikes && os.needed_lfp
            declone(Jparams_spike.f,Jparams_lfp.f);
            declone(Jparams_spike.confC,Jparams_lfp.confC);
            declone(Jparams_spike.N,Jparams_lfp.N);
            declone(Jparams_spike.nfft,Jparams_lfp.nfft);
        end
        T = [];
    else
        % Precalculates both time series and fft data
%         profile on
        [spike_all0,lfp_all0,T,J0_spike,J0_lfp,Jparams] = preload_and_save_tiled(os,md);
%         profile viewer
%         profile off
%         spike_all0 = [];
%         lfp_all0 = [];
%         J0_spike = [];
%         J0_lfp = [];
%         Jparams = [];
    end
    toc(ticID);
    
    % Pack everything into a single structure
    alldata.spike_all0 = spike_all0;
    alldata.lfp_all0 = lfp_all0;
    alldata.T = T;
    alldata.J0_spike = J0_spike;
    alldata.J0_lfp = J0_lfp;
    alldata.Jparams = Jparams;
    alldata.bad_trials = bad_trials;
    alldata.outname = outname;
    
    calc_stats(os,md,alldata);
end

end


function calc_stats(os,md,alldata)
    %% function calc_stats(filename_orig, sfc_mode, curr_stage, bad_trials, filenum, i0, j0, outname_embed_i0j0)
    
    
    %calc_stats(filename_orig, sfc_mode, curr_stage, bad_trials, filenum, i0, j0, outname_embed_i0j0)
%     profile on
    ticID2 = tic;
    % Pull out from os
    sfc_mode = os.sfc_mode;
    curr_stage = os.curr_stage;
    permutation_test = os.permutation_test;
    outname_embed_i0j0 = os.outname_embed_i0j0;
    filenum = os.filenum;
    filename_orig = os.fname;
    i0 = os.i0;
    outname = alldata.outname;
    
    % Extract SFC name stuff
    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix, do_adjacent, ctgsetli_mode, thinning_mode, tapers_mode, baseline_subtract, permutation_test, ue_pairs, coh_debias_mode, do_partial_coherence, ctgsetli_mode2, Nwind_mode] = build_sfcmode(sfc_mode, mode_subgroups);
    
    % Pull out from alldata
    bad_trials=alldata.bad_trials;
    
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
    
    % Calculate ctgsetli0
    [ctgsetli0,N_ctgs_base,N_ctgs_extras] = get_ctgsetli(sfc_mode,md,ctgsetli_mode,ctgsetli_mode2,filenum);
    N_categories = size(ctgsetli0,2);    % Update N_categories
    if strcmp(filename_orig,'L112106.mat') && (ue_pairs ~= 5) && (ue_pairs ~= 7)
        Ntrials = size(alldata.lfp_all0,2);
        ctgsetli0 = ctgsetli0(1:Ntrials,:);       % Shorten ctgsetli if we're missing lfp data. Should only be for file #29. Only do this if we are consider LFP
    end
    
    % Calculate ctgsetli maxes 
    ctgsetli_maxes = calc_ctgsetli_maxes(ctgsetli0,ctgsall,coh_debias_mode);

    % Assign function handle for spike-field or field-field coherence
    coherence_func = @calc_sfc_struct;
    
    % Do analysis of coherence! (permutation test or otherwise)
    if permutation_test == 0    % Default mode
        permutation_mode = 0;
        sfc = coherence_func(os,md,alldata, ctgsetli0, ctgsetli_maxes,i0);
        
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
            
            ctgsetli_curr = ctgsetli0;
            
            sfc_stat = coherence_func(os,md,alldata, ctgsetli_curr, ctgsetli_maxes,i0curr(1:2));
            isphi = isfield(sfc_stat{end},'phiave');
            if isphi; isphi = ~isempty(sfc_stat{end}.phiave); end
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
            Npermut = 150;
                if floor(sfc_mode) == 44
                    Npermut = 50;          % Do fewer permutations when in cfc mode for calculating PAC
                end
                if os.is_spectrogram
                    Npermut = 50;
                end
            % Npermut = 2;
            dCave_null = zeros([sz, Npermut]);
            Cavep1 = zeros([sz, Npermut]);
            Cavep2 = zeros([sz, Npermut]);
            phip1 = zeros([sz, Npermut]);
            phip2 = zeros([sz, Npermut]);
            permutation_mode = permutation_test;    % Turn on the appropriate permutation test (there is a handy 1:1 map between the test and the mode)
            clear sfc_temp_all sfc_temp
            for k = 1:Npermut
                ctgsetli_curr = permute_ctgsetli(ctgsetli0, i0curr(1:2));
                sfc_temp_all{k} = coherence_func(os,md,alldata, ctgsetli_curr, ctgsetli_maxes,i0curr(1:2));
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
            sfc{pi2}.Npermut = Npermut;
        end
        
        % Singleton fields
        singletonfieldnames = sfc_temp{1}.singletonfieldnames;
        for i = 1:length(singletonfieldnames)
            if isfield(sfc_temp{1},singletonfieldnames{i})
                sfc{1}.(singletonfieldnames{i}) = sfc_temp{1}.(singletonfieldnames{i});
            end
        end
        

        % Merge permutation tests!
        sfc = permute_null_to_pvals(sfc);
    end
    
    toc(ticID2);
%     profile viewer
%     profile off
    save(outname,'sfc');
   
end


function sfc = calc_sfc_struct(os,md,alldata, ctgsetli, ctgsetli_maxes,i0)
    %% sfc = calc_sfc_struct(lfp_sample, lfp_sample_indices, filename_orig, sfc_mode, curr_stage, bad_trials, ctgsall, filenum, i0, j0, permutation_mode)
    
    % Pull options
    % (alternatively, use vars_pull(os))
    filenum = os.filenum;
    filename_orig = os.fname;
    sfc_mode = os.sfc_mode;
    mode_subgroups = os.mode_subgroups;
    sfc_mode_group = os.sfc_mode_group;
    
    %i0 = os.i0;
    j0 = os.j0;

    exclude_clipping_trials = os.exclude_clipping_trials;
    curr_stage = os.curr_stage;
    Nwind = os.Nwind;
    dt = os.dt;
    thresh = os.thresh;
    params = os.params;

    
    % Pull alldata
    bad_trials = alldata.bad_trials;
    
    % Function specific options
    plot_debug1_thinning = 0;
    StarRate = get_rateStar;
    
    
    % Pull filename info (note - these variables are already pulled in get_options; I am just loading them again for simplicity
    [fname_suffix, do_adjacent, ctgsetli_mode, thinning_mode, tapers_mode, baseline_subtract, permutation_test, ue_pairs, coh_debias_mode, do_partial_coherence, ctgsetli_mode2, Nwind_mode] = build_sfcmode(sfc_mode, mode_subgroups);
    
    if sfc_mode_group == 31;    % Can't do baseline subtract, since we're using average signals here.
        error('Baseline subtract was turned on for correlation study \n');
    end
    
    
    N_categories = size(ctgsetli,2);
    Nunits = length(md.unit_names);    % Number of units tracked
    Nelects = length(md.lfp_names);      % Number of electrodes
    
    
    [mypairs, Ncoherences] = build_mypairs(ue_pairs,md);
    

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
        trials = find(ctgsetli(:,i));
        ctgsetli_max_curr = ctgsetli_maxes(i);
        
        ticID = tic;
        
        if i > N_categories; continue; end
        switch sfc_mode_group
            
            case {2, 3, 22, 23, 41, 45, 52}
                sfc{i}.Cave = [];
                sfc{i}.S12ave = [];
                sfc{i}.S1ave = [];
                sfc{i}.S2ave = [];
                sfc{i}.phiave = [];
                sfc{i}.confC = [];
                sfc{i}.phistd = [];
                sfc{i}.Cerr1 = [];
                sfc{i}.Cerr2 = [];
        end
        
        
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
        

        switch sfc_mode_group
            case {2,3,22,23,41,45}         % SFC / FFC, spectrogram or not
                
                
                func_data.os = os;
                func_data.alldata = alldata;
                func_data.mypairs = mypairs;
                [Cave,phi,S12ave,S1ave,S2ave,f,confC]=trials_downsample(ctgsetli_max_curr,trials,@postcalc_coherency,func_data);
                if any(sfc_mode_group == [41,45])
                    phi = [];   % Drop phi since it's useless to save space
                end

                sfc{i}.Cave = [sfc{i}.Cave Cave];
                sfc{i}.S12ave = [sfc{i}.S12ave S12ave];
                sfc{i}.S1ave = [sfc{i}.S1ave S1ave];
                sfc{i}.S2ave = [sfc{i}.S2ave S2ave];
                sfc{i}.phiave = [sfc{i}.phiave phi];
                sfc{1}.f = f;
                sfc{1}.f2 = alldata.T;
                sfc{i}.confC = [sfc{i}.confC];
                sfc{i}.phistd = [];

                sfc{i}.Cerr1 = [];
                sfc{i}.Cerr2 = [];
                

                
            case 52
                    spike_all0 = alldata.spike_all0(:,trials,:);
                    switch tapers_mode
                        case 0          % Average entire spike train
                            
                            % % More complicated way of doing things % %
                            % spikerates2 = currspike1_ctg / dt; % Converts spike event (0 or 1) into spike rate per time bin.
                            % spikerates2 = mean(spikerates2,2); % Mean spike rate across trials
                            % spikerates2 = mean(spikerates2(:)); % Mean spike rate across all time.
                            
                            % % % This should be equal to the above % % %
                            spikerates2 = squeeze(mean(sum(spike_all0,1) / (size(spike_all0,1)*dt))); % Should work out to the same as this.
                            sfc{i}.Cave = [sfc{i}.Cave repmat(spikerates2',5,1)]; 
                            sfc{1}.f = 1:5;
                        case {1,2,3,4,5}
                            Nwind = os.Nwind;
                            spikerates2 = squeeze(mean(spike_all0,2)) / dt;
                          
                            spikerates2 = sgolayfilt(spikerates2,3,Nwind);   % Apply 3rd-order filter, 151 ms 

                            sfc{i}.Cave = [sfc{i}.Cave spikerates2];
                            sfc{1}.f = [1:length(spikerates2)]*dt;
                    end
        end

        % Add some i0 specific metadata
        sfc{i}.Ntraces = [sfc{i}.Ntraces length(trials)];
        sfc{i}.mypairs = [sfc{i}.mypairs mypairs(:,:)'];
        
        sfc{i}.sum_ctgsetli_maxes = [sfc{i}.sum_ctgsetli_maxes ctgsetli_max_curr];
        
        

        if any(ue_pairs == [2,3,5,7])           % Only save spike info if we're not in LFP only mode
            
            if os.is_spectrogram; spike_all0_1 = alldata.spike_all0(:,trials,:,mypairs(:,1));
            else spike_all0_1 = alldata.spike_all0(:,trials,mypairs(:,1));
            end
            sfc{i}.Ntraces_below_thresh = [sfc{i}.Ntraces_below_thresh squeeze(sum(     mean(spike_all0_1,1)    /dt < thresh, 2))];
            sfc{i}.spikerate_mu = [sfc{i}.spikerate_mu squeeze(mean(mean(spike_all0_1,1),2))/dt];
            
            if ue_pairs == 5
                if os.is_spectrogram; spike_all0_2 = alldata.spike_all0(:,trials,:,mypairs(:,2));
                else spike_all0_2 = alldata.spike_all0(:,trials,mypairs(:,2));
                end
                sfc{i}.Ntraces_below_thresh2 = [sfc{i}.Ntraces_below_thresh2 squeeze(sum(     mean(spike_all0_2,1)    /dt < thresh, 2))];
                sfc{i}.spikerate_mu2 = [sfc{i}.spikerate_mu2 squeeze(mean(mean(spike_all0_2,1),2))/dt];    % total # spikes / total time
            end
        end

        

        sfc{i}.elapsedTime = [sfc{i}.elapsedTime, toc(ticID)];


    end
    
    
    % Save output, unless we're in raster plotting mode.
    sfc{1}.singletonfieldnames = {'singletonfieldnames','f','f2','f3' ...
        'F','F2','T', ...                                                   % For 3D plots (coherograms)
        'unit_names','filename','sum_ctgsetli',...
        'include_switch_trials','N_ctgs_base','N_ctgs_extras','stagesir', ...
        'md_lfp_names','md_unit_names',...
        'spikethresh','StarRate','params','i0','j0','os'};
    sfc{1}.unit_names = unit_names;
    sfc{1}.filename = filename_orig;
    sfc{1}.sum_ctgsetli = sum(ctgsetli);
    sfc{1}.include_switch_trials = ctgsetli_mode; % Metadata
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
    sfc{1}.os = os;
    


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
    
    % Make sure the two original groups are non-overlapping
    N_overlapping = sum(ctgs1 & ctgs2);
    if N_overlapping > 0
        error('Error, cannot shuffle groups. Groups contain overlapping items.');
    end
    
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


function [mypairs,Ncoherences] = build_mypairs(ue_pairs,md)
    %% function mypairs = build_mypairs(ue_pairs,md)
    % Set up mypairs for iteration

    Nunits = length(md.unit_names);    % Number of units tracked
    Nelects = length(md.lfp_names);      % Number of electrodes
    
    switch ue_pairs
        case {0,2};                                             % Current unit versus current/adjacent field (SFC)
            Ncoherences = Nunits;
            
            units1 = [1:Nunits];
            
            if ue_pairs == 0;
                do_adjacent = 0;
            elseif ue_pairs == 2;
                do_adjacent = 1;
            end
            
            lfp2 = get_corresponding_electrodes_vectorized(units1,md,do_adjacent);
            
            mypairs = [units1; lfp2]';
            clear units1 lfp2
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
end


function temporary_code
  %% TEmp code
  
    
    if permutation_mode == 1
        if length(i0) ~= 2; warning('When permutation_mode=1, length(i0) must be 2.'); return; end
        ctgsetli_good{j} = permute_ctgsetli(ctgsetli_good{j}, i0);
    end
        

end


