

%% Set up everything


format compact
format long g

run_setdefaultfig

%     load paths
    addpath(genpath('./funcs_supporting_local'));
    addpath(genpath('./funcs_run_analysis_precalc'));

    
    sfc_mode =  45.6013111011; 
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


    file_range = 1:Nfiles;
%    file_range = [4];
%    file_range = [44];
%    file_range = [70];
%     file_range = [1];
%     file_range = [1:Nfiles];

filenum = 1; 

% File name
fname = file_list{filenum};
fprintf ('Processing file %d of %d, name %s \n',filenum,file_range(end),fname);

% Load options structure
os = get_options(sfc_mode, curr_stage,i0,j0,outname_embed_i0j0,file_list,filenum);

%%

[spike_all0,lfp_all0,T,J0_spike,J0_lfp,Jparams] = preload_and_save_tiled(os,md);



%% Load data

md = load('/Users/davestanley/Other/data/roy/ctgX_merged_final/L011107.mat');

load('/Users/davestanley/Other/data/roy/CRC_Miller_data/intmatrix/L011107.mat')

lfpmat = double(lfpmat_int16); clear lfpmat_int16


%%



%% 

[~, mode_subgroups] = decode_sfc_mode(sfc_mode);
[fname_suffix, do_adjacent, ctgsetli_mode, thinning_mode, tapers_mode, baseline_subtract, permutation_test, ue_pairs, coh_debias_mode, do_partial_coherence, ctgsetli_mode2, Nwind_mode] = build_sfcmode(sfc_mode, mode_subgroups);
    
[ctgsetli0,N_ctgs_base,N_ctgs_extras] = get_ctgsetli(sfc_mode,md,ctgsetli_mode,ctgsetli_mode2,1);


%%

vars_pull(os);

[Ssp, phi, S12ave, S1ave, S2ave, F, confC, phistd, Cerr, T] = func2spectrogram(currlfp1_ctg,currlfp2_ctg,'fs',1/dt,'Nwind',Nwind,'fract_overlap',fract_overlap,'fname',fname);



