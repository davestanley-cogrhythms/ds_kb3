function run_analysis_FR (sfc_mode, stage, file_range)

% % % Test Fries method on Kyle's simulated data. Then test Kyle's method
% and compare
% % % 
format compact
format long g

run_setdefaultfig
addpath(genpath('./funcs_supporting_local'));

if ~exist('sfc_mode','var'); sfc_mode = 5.0; end
    % sfc_mode = 5 -> unit analysis             % MUST BE 5 for this code!
        % sfc_mode = 5.1 -> unit analysis, only this version reports the smoothed mean firing rates instad of firing rates per trial
    
    
if ~exist('stage','var'); stage = 2; end
    % stage = 1 - pre sample on
    % stage = 2 - sample stage
    % stage = 3 - delay stage
    % stage = 4 - sample + delay stage
    % stage = 5 - all
    % stage = 6 - Test1 - match or non-match
    % stage = 7 - TDelay - test delay
    % stage = 8 - Test2 - match (if prior was non-match)
    


% path_metadata = getpath('path_metadata');
metadata_path = fullfile(getpath('path_metadata'));
% path_lfp = getpath('path_lfp');


file_list = get_filelist(metadata_path);
Nfiles = length(file_list);

% Load badtrials
bad_trials = prep_badtrials(stage);

if  ~exist('file_range','var')
    file_range = 1:Nfiles;
%     file_range = [25:30];
    file_range = [1:Nfiles];
end

for i = file_range
    fprintf ('Processing file %d of %d, name %s \n',i,file_range(end),file_list{i});
    calc_stats(file_list{i},sfc_mode,stage,i);
end

end



function calc_stats(filename, sfc_mode, curr_stage, filenum)

    
    plot_on = 0;
    plot_debug = 0;
    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix] = build_sfcmode(sfc_mode, mode_subgroups);
    sfc_mode_group = floor(sfc_mode); % Mode group might be more than 1 dimension, so use floor instead (i.e. sfc_mode = 23.2 would yield mode_group == [2 3] instad of 23)

    if sfc_mode_group == 5 && mode_subgroups(1) == 1    % Recall; mode=5 and mode_subgroupps=1 corresponds to plotting smoothed spk traces
        curr_stage = 5; % For looking at whole set of activity
        apply_filter = 1;   % Apply sgolay filter to spiking data
    end
    
    dt = get_dt;
    
    outpath = fullfile(getpath('path_buffer_curr'),['mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(curr_stage)]);
    outname = fullfile(outpath,filename);
    mkdir_dave(outpath);
    if exist(outname,'file')
        fprintf('File %s exists...skipping\n',outname);
        return
    end

    md = load (fullfile(getpath('path_metadata'),filename));     % Load metadata
    
    Ntrials = length(md.start_times);  % Total number of trials
    Nsamples = length(md.sample_on_times); %Number of "good" samples - when trial reached sample_on stage
    Nunits = length(md.unit_names);    % Number of units tracked
    
    ctgsetli_mode = 0;
    ctgsetli_mode2 = 0;
    warning('This is a hack! ... but this code is outdated anways.');
    [ctgsetli,N_ctgs_base,N_ctgs_extras] = get_ctgsetli(sfc_mode,md,ctgsetli_mode,ctgsetli_mode2, filenum);
    
    N_categories = size(ctgsetli,2);    % Update N_categories
        

    cat_range = 1:(N_categories);          % 9th category is all data!
    
    for i = cat_range
        switch sfc_mode_group
            case 5
            sfc{i}.spikerates = [];
        end

        %for j = [18]
        for j = 1:Nunits
            curr_unit = j;
            unit_name = md.unit_names{curr_unit}; unit_name= strrep(unit_name,'_',' ');

           
            switch mode_subgroups(1)    % Check if mode is to do unit analysis trace plotting
                case 0
                    spikerates = get_numspikes_fromfile(curr_unit,curr_stage,md);
%                         % Normalizing
%                         numspikes = (numspikes - min(numspikes)) / (max(numspikes) - min(numspikes));
                    spikerates = spikerates(ctgsetli(:,i));
                case 1
                    currspike = get_spike_ts(curr_unit,curr_stage,md);
                    currspike = currspike(:,ctgsetli(:,i));
                    spikerates = mean(currspike,2); spikerates = spikerates / dt;
                    if apply_filter
                        spikerates = sgolayfilt(spikerates,3,round(0.151/dt));   % Apply 3rd-order filter, 151 ms 
                    end
            end

            
            
            switch sfc_mode_group
                case 5
                sfc{i}.spikerates = [sfc{i}.spikerates spikerates(:)];
            end
            
        end
    end
    
    % Save output, unless we're in raster plotting mode.
    sfc{1}.unit_names = md.unit_names;
    sfc{1}.filename = filename;
    sfc{1}.tapers = get_tapers;
    sfc{1}.sum_ctgsetli = sum(ctgsetli);
    sfc{1}.N_ctgs_base = N_ctgs_base;
    sfc{1}.N_ctgs_extras = N_ctgs_extras;
    sfc{1}.stagesir = get_stagesir(curr_stage);
    save(outname,'sfc');
    
end


