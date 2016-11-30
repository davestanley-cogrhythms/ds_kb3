
%% Test2
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
if ~exist('sfc_mode','var'); sfc_mode = 45.6013111011; end
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
    
if ~exist('curr_stage','var'); curr_stage = 2; end
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
if ~exist('i0','var'); i0 = [1:12]; end       % There are a max of 11 categories Roy data; There are fewer for Cromer
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

for filenum = 1
    
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
