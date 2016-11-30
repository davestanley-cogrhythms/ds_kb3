machine='dave';

% % Get home folder
if ispc; userdir= getenv('USERPROFILE'); 
else userdir= getenv('HOME'); 
end

% % System-specific customization


% Is this machine the one you sold? -Kelly
% % Dave's Thinkpad T520 configurations (comment / uncomment)
% path_src = fullfile(filesep,'data','Anco','src');
% path_roy = fullfile(filesep,'media','davestanley','Dave','Users','Computer','data','roy');

switch machine
  case 'dave'
    % Dave's SCC configurations (comment / uncomment)
    path_src = fullfile(userdir,'src'); % Path containing all code
    %path_roy = fullfile(userdir,'Miller','roy_preplexon'); % Path containing all data prior to plexon filtering
    path_roy = fullfile(userdir,'Miller','roy'); % PAth containing all data
    path_synth = fullfile(userdir,'Miller','synth');    % Path for synthetic data
  case 'kelly'
    % Kelly's SCC configurations (comment / uncomment)
    path_src = fullfile(userdir,'src'); % Path containing all code
    path_roy = fullfile(userdir,'crc-nak','roy'); % Path containing all data
  otherwise
    disp(['The environment ' machine ' is not defined.']);
end

% % SHOULD NOT HAVE TO CHANGE BELOW HERE
% % Paths from git repositories (should not require editing)
basepath = fullfile(path_src,'ds_kb');
chronuxpath = fullfile(path_src,'chronux');
sfcpath = fullfile(path_src,'spike_field_assoc_dev');

% Paths to data and metadata (should not require editing)
path_roydata_orig = fullfile(path_roy,'CRC_Miller_data');
path_metadata = fullfile(path_roy,'ctgX_merged_final');
path_conditions = fullfile(path_roy,'conditions_codes');  % Path to data that maps condition codes onto test image ctgs
path_lfp = fullfile(path_roydata_orig,'intmatrix');
path_lfp_sample_delay = fullfile(path_roy,'lfp_sample_delay');
path_preferred = fullfile(path_roy,'new_preferred_roy');
path_sample_on_times_per_trial = fullfile(path_roy,'sample_on_times_per_trial'); % For Kelly's samples code
path_electrode_locations = fullfile(path_roy,'electrode_locations');

% Output paths
path_SFC = fullfile(path_roy,'SFC_out');        % I am no longer using this... will remove soon - Kelly, do you use this path?

% Paths to elsewhere
path_roy_newctgX = fullfile(path_roy,'new_ctgX_files');


% Buffer for run_analysis.m output
path_buffer = fullfile(path_roy,'buffer');                              % Outputs of run_analysis
path_buffer_curr = fullfile(path_buffer,'analyze_data7b');              

% Cache for run_analysis_precalc.m
path_precalcD = fullfile(path_roy,'precalcD');
path_precalcD_curr = fullfile(path_precalcD,'analyze_data7b');
path_precalcJ = fullfile(path_roy,'precalcJ');
path_precalcJ_curr = fullfile(path_precalcJ,'analyze_data7b');


% Buffer for running plots_preferred
path_buffer_PL = fullfile(path_roy,'buffPL');                           % Things used in plots_preferred (matifles buffer, badchannels, etc)
path_buffer_currPL = fullfile(path_buffer_PL,'analyze_data7');            
path_bufffer_SFCsurvey = fullfile(path_buffer_currPL,'betaSFC_survey');
path_buffer_mypreferred = fullfile(path_buffer_currPL,'preferred_units');
path_buffer_specialization_stats = fullfile(path_buffer_currPL,'stats_specialize');
path_buffer_PEV = fullfile(path_buffer_currPL,'stats_PEV');
path_matfiles_store = fullfile(path_buffer_currPL,'matfiles_store');
path_badchannels = fullfile(path_buffer_currPL,'badchannels_list');

% Analysis data
% path_run_analysis2 = fullfile(path_buffer,'run_analysis2');

% Path for synthetic data
path_synth_buffer = fullfile(path_synth,'buffer1');
path_synth_synthout = fullfile(path_synth_buffer,'synth_out');
path_synth_matfiles = fullfile(path_synth_buffer,'matfiles');



