

function run_analysis_survey_spectrogram_vs_ev (sfc_mode, stage, file_range)

% % % Test Fries method on Kyle's simulated data. Then test Kyle's method
% and compare
% % % 
format compact
format long g

run_setdefaultfig
addpath(genpath('./funcs_supporting_local'));

if ~exist('sfc_mode','var'); sfc_mode = 8.0; end
    % sfc_mode = 1 -> fries
    % sfc_mode = 2 -> chronux
    % sfc_mode = 3 -> cohernece_spect_wrapper (segmented)
    % sfc_mode = 4 -> SFC analysis trial-by-trial
    % sfc_mode = 5 -> unit analysis
        % sfc_mode = 5.1 -> unit analysis, only this version reports the smoothed mean firing rates instad of firing rates per trial
    % sfc_mode = 6 -> raster_plotting
    % sfc_mode = 7 -> Lepage (not yet implemented)
    % sfc_mode = 8 -> spectrogram
    % sfc_mode = 9 -> evoked responses
    
    % sfc_mode = *.2   -> use instad adjacent electrode! - not applicable for unit analysis (sfc_mode=5.*)
    % sfc_mode = *.?1  -> include analysis of only switch trials. Doubles the size of output, usually.
if ~exist('stage','var'); stage = 3; end
    % stage = 1 - pre sample on
    % stage = 2 - sample stage
    % stage = 3 - delay stage
    % stage = 4 - sample + delay stage
    % stage = 5 - all


% path_metadata = getpath('path_metadata');
path_lfp_sample_delay = getpath('path_lfp_sample_delay');
% path_lfp = getpath('path_lfp');


file_list = get_filelist(path_lfp_sample_delay);
Nfiles = length(file_list);

% Load badtrials
bad_trials = prep_badtrials(stage);

if  ~exist('file_range','var')
    file_range = 1:Nfiles;
%     file_range = [25:30];
    file_range = [2];
end

for i = file_range
    fprintf ('Processing file %d of %d, name %s \n',i,length(file_list),file_list{i});
    calc_stats(file_list{i},sfc_mode,stage,bad_trials,i);
end

end



function calc_stats(filename, sfc_mode, curr_stage, bad_trials, filenum)

    if ~exist('filename','var'); filename = 'L091906.mat'; end
    if ~exist('sfc_mode','var'); sfc_mode = 2; end
                % sfc_mode = 1 -> fries
                % sfc_mode = 2 -> chronux
                % sfc_mode = 3 -> Lepage (not yet implemented)
                % sfc_mode = 4 -> mtcoherencysegc (segmented)

    % Note: indices refers to indices of lfp matrix
    % Trials refers to indices of the ctgX_trials arrays
    
    plot_on = 0;
    plot_debug = 0;
    [mode_group, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix, include_switch_trials, do_adjacent] = build_sfcmode(sfc_mode, mode_subgroups);
    sfc_mode_group = floor(sfc_mode); % Mode group might be more than 1 dimension, so use floor instead (i.e. sfc_mode = 23.2 would yield mode_group == [2 3] instad of 23)
    

    exclude_clipping_trials=1;
    loadall_lfps = 1;      % Used for comparing the current spiking sequence against all LFP's              
    
    if sfc_mode_group == 3 || sfc_mode_group == 8 || sfc_mode_group == 9     % For spectrograms and evoked responses, may as well look at all activity
        curr_stage = 5;
        Nwind = round(0.5/get_dt);
    end
    
    if sfc_mode_group == 6     % Raster plotting mode
        curr_stage = 5; % For looking at the whole set of activity.
        plot_sfc_on = 1;
    end
    
    if sfc_mode_group == 5 && mode_subgroups(1) == 1    % Recall; mode=5 and mode_subgroupps=1 corresponds to plotting smoothed spk traces
        curr_stage = 5; % For looking at whole set of activity
        apply_filter = 1;   % Apply sgolay filter to spiking data
    end
    
    dt = get_dt;
    wind_s = round(0.2/dt);     % Window size for Fries method
    thresh = 5;                 % Firing rate threshold for rejecting traces
    
    
    outpath = fullfile(getpath('path_buffer_curr'),['mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(curr_stage)]);
    outname = fullfile(outpath,filename);
    if sfc_mode_group ~= 6; mkdir_dave(outpath); end
    if exist(outname,'file')
        fprintf('File %s exists...skipping\n',outname);
%         return
    end

    md = load (fullfile(getpath('path_metadata'),filename));     % Load metadata
    %load (fullfile(getpath('path_lfp'),filename));     % Load full file
    %lfpmat = double(lfpmat_int16); clear lfpmat_int16;
    
    if sfc_mode_group ~= 5
        load (fullfile(getpath('path_lfp_sample_delay'),filename));
        lfp_sample = double(lfp_sample);
        %lfp_indices = recall_lfpir(lfp_sample_indices); % I don't think I use this
        lfp_indices_abs = bounds_to_indices([lfp_sample_indices(:,1) lfp_sample_indices(:,3)]');
        %lfp_indices_abs2 = recall_lfpia(lfp_sample_indices);
    end
    
    
    
    Ntrials = length(md.start_times);  % Total number of trials
    Nsamples = length(md.sample_on_times); %Number of "good" samples - when trial reached sample_on stage
    Nunits = length(md.unit_names);    % Number of units tracked
    Nelects = length(md.lfp_names);      % Number of electrodes
    
    
    if include_switch_trials
        
        tsw = get_switch_trials(md,Ntrials);
        
    end
    
    
    N_categories = 4*2 + 1;
    ctgsetli = false(length(get_good_samples(md.ctg1_trials,md)),N_categories);
    ctgsetli(:,1) = get_good_samples(md.ctg1_trials,md);
    ctgsetli(:,2) = get_good_samples(md.ctg2_trials,md);
    ctgsetli(:,3) = get_good_samples(md.ctg3_trials,md);
    ctgsetli(:,4) = get_good_samples(md.ctg4_trials,md);
    
    ctgsetli(:,5) = get_good_samples(md.ctg1_nr_trials,md);
    ctgsetli(:,6) = get_good_samples(md.ctg2_nr_trials,md);
    ctgsetli(:,7) = get_good_samples(md.ctg3_nr_trials,md);
    ctgsetli(:,8) = get_good_samples(md.ctg4_nr_trials,md);
    ctgsetli(:,9) = ones(size(ctgsetli,1),1);
    
    if include_switch_trials
        tsw_mat = repmat(tsw(:),1,size(ctgsetli,2));
        ctgsetli = [ctgsetli [ctgsetli & tsw_mat]];   % Double size of ctgsetli
    end
    N_categories = size(ctgsetli,2);
    
    
    %[in_ctg in_ctgnr] =  ctg_query(2,md);
    


    % Plots the default time series with the "samples" overlaid
    if plot_debug
%         curr_unit = 1;
%         [currlfp currspike lfpia_stage] = get_lfp_and_spike_pairs(lfp_sample,curr_unit,curr_stage,md,lfp_sample_indices);
        %plot_abs_vs_relative_lfp_and_spikes(filename,lfp_sample_indices,curr_unit,md,currlfp,currspike,lfpia_stage);
    end
    
    
    
    i=1; j=1;
    C_arr = [];
    %pause;
    if sfc_mode_group == 6
        cat_range = 9;                          % Only look at the 9th category!
    else
        cat_range = 1:(N_categories);          % 9th category is all data!
        %cat_range = 1;
    end
    
    for i = cat_range
        switch sfc_mode_group
            case 1
                sfc{i}.Cave = [];       % Make it a cell array to allow for variation amongst # of categories
                sfc{i}.lfp_wind = [];       % Would be useful if we ever want to store individual windows

            case 2
                sfc{i}.Cave = [];
                sfc{i}.S12ave = [];
                sfc{i}.S1ave = [];
                sfc{i}.S2ave = [];
                sfc{i}.phiave = [];
            
            case 3
                sfc{i}.Cave = [];
                sfc{i}.S12ave = [];
                sfc{i}.S1ave = [];
                sfc{i}.S2ave = [];
                
            case 4
                sfc{i}.S12ave = [];
                sfc{i}.S1ave = [];
                sfc{i}.S2ave = [];

            case 5
                sfc{i}.numspikes = [];
                
            case 8
                sfc{i}.S = [];
            case 9
                sfc{i}.E = [];
        end
%         for j = 1:Nunits
        for j = [8]
            i;
            curr_unit = j;
            unit_name = md.unit_names{curr_unit}; unit_name= strrep(unit_name,'_',' ');
            adj_missing = 0;

            if sfc_mode_group ~=5
                %[currlfp_old currspike_old lfpia_stage_old] = get_lfp_and_spike_pairs(lfp_sample,curr_unit,curr_stage,md,lfp_sample_indices);  % Old code - I've split this up below
                [currlfp, lfpia_stage, adj_missing] = get_unitLFP_ts(lfp_sample,curr_unit,curr_stage,md,lfp_sample_indices,do_adjacent);
                if loadall_lfps; [currlfp_all,lfpia_stage] = get_LFP_ts(lfp_sample,curr_stage,md,lfp_sample_indices); end
                [currspike] = get_spike_ts(curr_unit,curr_stage,md);
                if plot_debug
                    plot_abs_vs_relative_lfp_and_spikes(filename,lfp_sample_indices,curr_unit,md,currlfp,currspike,lfpia_stage);
                end
                
                if strcmp(filename,'L112106.mat')
                    ctgsetli = ctgsetli(1:size(currlfp,2),:);       % Shorten ctgsetli if we're missing lfp data. Should only be for file #29.
                end
                
                if exclude_clipping_trials
                    curr_badtrials = get_curr_badtrials(bad_trials,filenum,curr_unit,md);
                else
                    curr_badtrials = false(Nsamples,1); % All trials are good!
                end
                
                currlfp_ctg = currlfp(:,ctgsetli(:,i) & ~curr_badtrials);
                currspike_ctg = currspike(:,ctgsetli(:,i) & ~curr_badtrials);
                if loadall_lfps; currlfp2_ctg = currlfp_all(:,ctgsetli(:,i) & ~curr_badtrials,:); end
                
                
            else
                switch mode_subgroups(1)    % Check if mode is to do unit analysis trace plotting
                    case 0
                        numspikes = get_numspikes_fromfile(curr_unit,curr_stage,md);
%                         % Normalizing
%                         numspikes = (numspikes - min(numspikes)) / (max(numspikes) - min(numspikes));
                        numspikes = numspikes(ctgsetli(:,i));
                    case 1
                        [currspike] = get_spike_ts(curr_unit,curr_stage,md);
                        currspike = currspike(:,ctgsetli(:,i));
                        numspikes = mean(currspike,2); numspikes = numspikes / dt;
                        if apply_filter
                            numspikes = sgolayfilt(numspikes,3,round(0.151/dt));   % Apply 3rd-order filter, 151 ms 
                        end
                end
                        
                    
            end
            
            
            switch sfc_mode_group
                case 1
                    %[C ftemp] = sfc_fries2D(currspike(:,1),currlfp(:,1),dt,wind_s);
                    
                    %tic
                    [Cave f lfp_wind] = sfc_fries2D(currspike_ctg,currlfp_ctg,dt,wind_s);
                    lfp_wind_mu = mean(lfp_wind,2);
                    lfp_wind_std = std(lfp_wind,[],2);
                    lfp_wind_ste = lfp_wind_std/sqrt(size(lfp_wind,2));
                    %toc
                    
                    sfc{i}.Cave = [sfc{i}.Cave Cave(:)];
                    sfc{i}.lfp_wind = [sfc{i}.lfp_wind lfp_wind_mu(:)];
                    sfc{i}.Nunits(j) = size(lfp_wind,2);
                    sfc{i}.spikerate_mu(j) = sum(sum(currspike_ctg,1)) / (size(currspike_ctg,1)*get_dt*sfc{i}.Ntraces(j));    % total # spikes / total time
                    %sfc{i}.spikes_per_trial{j} = sum(currspike_ctg,1)';
                    sfc{1}.f = f;
                    
                case 2
                    %tic
                    params.Fs = 1/dt;
                    params.tapers=get_tapers;

                    params.trialave = 1;
                    sfc{i}.Ntraces(j) = size(currspike_ctg,2);
                    sfc{i}.Ntraces_below_thresh(j) = sum(sum(currspike_ctg,1) < thresh);
                    sfc{i}.spikerate_mu(j) = sum(sum(currspike_ctg,1)) / (size(currspike_ctg,1)*get_dt*sfc{i}.Ntraces(j));    % total # spikes / total time
                    sfc{1}.spikethresh = thresh;
                    %sfc{i}.spikes_per_trial{j} = sum(currspike_ctg,1)';
                    
                    %currlfp = zscore(currlfp); % Note- I found that taking the zscore doesn't make any difference. I guess this is because things are normalized already by the SFC operation.
                    
                    % Remove traces below a set number of spikes
%                     [currlfp_ctg,currspike_ctg] = trim_lowspike_ctgs(currlfp_ctg,currspike_ctg,thresh);
                    [Cave,phi,S12ave,S1ave,S2ave,f]=coherencycpb(currlfp_ctg,currspike_ctg,params);
                    %toc
%                     figure('Position',[  1          26        1280         769]);
%                     range = 1
%                     subplot(221); title([filename ' unit # ' num2str(j) ' ' unit_name])
%                     hold on; plot(f,Cave(:,range));legend('Cave'); xlabel('f (Hz)'); ylabel('SFC');
%                     subplot(222); hold on; plot(f,S1ave(:,range));legend('S1ave'); xlabel('f (Hz)'); ylabel('Power');
%                     subplot(223); hold on; plot(f,S2ave(:,range));legend('S2ave'); xlabel('f (Hz)'); ylabel('Power');
%                     subplot(224); hold on; plot(sum(currspike_ctg)); legend('Tot spikes in each trial'); xlabel('Trial #'); ylabel('# Spikes');

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

                    
                    sfc{i}.Cave = [sfc{i}.Cave Cave(:)];
                    sfc{i}.S12ave = [sfc{i}.S12ave S12ave(:)];
                    sfc{i}.S1ave = [sfc{i}.S1ave S1ave(:)];
                    sfc{i}.S2ave = [sfc{i}.S2ave S2ave(:)];
                    sfc{i}.phiave = [sfc{i}.phiave phi(:)];
                    sfc{1}.f = f;
                    
                    
                    
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
                    params.tapers=get_tapers;

                    params.trialave = 1;
                    sfc{i}.Ntraces(j) = size(currspike_ctg,2);
                    sfc{i}.Ntraces_below_thresh(j) = sum(sum(currspike_ctg,1) < thresh);
                    sfc{i}.spikerate_mu(j) = sum(sum(currspike_ctg,1)) / (size(currspike_ctg,1)*get_dt*sfc{i}.Ntraces(j));    % total # spikes / total time
                    sfc{1}.spikethresh = thresh;
                    %sfc{i}.spikes_per_trial{j} = sum(currspike_ctg,1)';

                    % Remove traces below a set number of spikes
%                     [currlfp_ctg,currspike_ctg] = trim_lowspike_ctgs(currlfp_ctg,currspike_ctg,thresh);
                    [Cave, phi, S12ave, S1ave, S2ave, f, T] = coherence_spect_wrapper(currlfp_ctg,currspike_ctg,'fs',1/dt,'params',params,'Nwind',Nwind);
                    
                    % Trim Cave to only go up to 250 Hz
                    index = find(f <= 250);
                    f = f(index);
                    Cave = Cave(index,:);
                    S1ave = S1ave(index,:);

                    % Pack variables
                    Cave = permute(Cave,[1,3,2]);
                    S1ave = permute(S1ave,[1,3,2]);
                    sfc{i}.Cave = cat(2,sfc{i}.Cave,Cave);
                    %sfc{i}.S1ave = cat(2,sfc{i}.S1ave,S1ave);    % Don't save this stuff to save space

                    sfc{1}.f = f;
                    sfc{1}.T = T;
                    
                case 4

                    params.Fs = 1/dt;
                    params.tapers=get_tapers;

                    params.trialave = 0;
                    
                    sfc{i}.Ntraces(j) = size(currspike_ctg,2);
                    sfc{i}.Ntraces_below_thresh(j) = sum(sum(currspike_ctg,1) < thresh);
                    sfc{i}.spikerate_mu(j) = sum(sum(currspike_ctg,1)) / (size(currspike_ctg,1)*get_dt*sfc{i}.Ntraces(j));    % total # spikes / total time
                    sfc{1}.spikethresh = thresh;
                    %sfc{i}.spikes_per_trial{j} = sum(currspike_ctg,1)';
                    
                    %currlfp = zscore(currlfp); % Note- I found that taking the zscore doesn't make any difference. I guess this is because things are normalized already by the SFC operation.
                    
                    % Remove traces below a set number of spikes
%                     [currlfp_ctg,currspike_ctg] = trim_lowspike_ctgs(currlfp_ctg,currspike_ctg,thresh);
                    [Cave,phi,S12ave,S1ave,S2ave,f]=coherencycpb(currlfp_ctg,currspike_ctg,params);
                    
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
                    sfc{i}.numspikes = [sfc{i}.numspikes numspikes(:)];
                    
                case 6
                    plot_unit_summary(currspike,curr_stage,ctgsetli,filename,md.unit_names{j},curr_unit,plot_sfc_on);
                    %plot_raster_with_sfc(currspike,curr_stage,ctgsetli,filename,md.unit_names{j});
                    
                case 8      % LFP Spectrogram
                    params.Fs = 1/dt;
                    params.tapers=get_tapers;

                    params.trialave = 1;
                    sfc{i}.Ntraces(j) = size(currspike_ctg,2);
                    [S, f, T] = spect_wrapper(currlfp_ctg,'fs',1/dt,'Nwind',Nwind,'params',params);
                    
                    % Trim Cave to only go up to 250 Hz
                    index = find(f <= 250);
                    f = f(index);
                    S = S(index,:);

                    % Pack variables
                    S = permute(S,[1,3,2]);
                    sfc{i}.S = cat(2,sfc{i}.S,S);

                    sfc{1}.f = f;
                    sfc{1}.T = T;

                    figure; plott_spect(currlfp_ctg,'fs',1/dt,'zscoreplot',1); ylim([0 250])
                    figure
                    mytrials = randperm(size(currlfp_ctg,2));
                    for iii = 1:size(currlfp_ctg,2)
                        clf;
                        currlfp_trial = currlfp(:,mytrials(iii));
                        subplot(211);
                        plott_fs(currlfp_trial,'fs',1/dt);xlabel('time(s)');
                        subplot(212); plott_spect(currlfp_trial,'fs',1/dt,'zscoreplot',1,'Nwind',round(0.5/dt)); ylim([0 250])
                        pause(0.2)
                    end
                    
                case 9      % Evoked potential
                    E = mean(currlfp_ctg,2);

                    % Pack variables
                    sfc{i}.E = cat(2,sfc{i}.E,E(:));

            end
            
            % Other stuff to save, independent of sfc_mode.
            sfc{1}.adj_was_missing(j) = adj_missing;
        end

    end
    
    if sfc_mode_group ~= 6
        % Save output, unless we're in raster plotting mode.
        sfc{1}.unit_names = md.unit_names;
        sfc{1}.filename = filename;
        sfc{1}.tapers = get_tapers;
        sfc{1}.ctgsetli = sum(ctgsetli);
        sfc{1}.do_adjacent = do_adjacent;   % Metadata
        sfc{1}.include_switch_trials = include_switch_trials; % Metadata
%         save(outname,'sfc');
    end
    
%     [pathstr, name, ext] = fileparts(filename);
%     %mkdir_dave(getpath('path_SFC'));
%     
%     switch sfc_mode_group
%         case 1
%             save(fullfile(getpath('path_SFC'),[name '_sfc_mode_group_' num2str(sfc_mode_group) '.mat']),'C_arr','f');
%         case 2
%             save(fullfile(getpath('path_SFC'),[name '_sfc_mode_group_' num2str(sfc_mode_group) '.mat']),'C_arr','S1_arr','f');
%             
%     end

end






function plot_abs_vs_relative_lfp_and_spikes(filename,lfp_sample_indices,curr_unit,md,currlfp,currspike,lfpia_stage)
        

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
%     params.tapers=[3 5];
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
    index = sum(currspike_ctg) >= thresh;
    currlfp_ctg = currlfp_ctg(:,index);
    currspike_ctg = currspike_ctg(:,index);
end




function curr_badtrials = get_curr_badtrials(bad_trials,filenum,curr_unit,md)
    curr_lfp = md.unit_LFP_pairs(2,curr_unit);
    curr_badtrials = bad_trials{filenum}(:,curr_lfp);
end



function tsw = get_switch_trials(md,Ntrials)
        
        for iii = 1:Ntrials
            [sch(iii)] = condition2morph(md.each_trial_condition(iii)); % Slow, but not limiting.
        end
        
        sch = (sch == 'A');
        sch_curr = sch(2:end);
        sch_prev = sch(1:end-1);

        tsw = sch_curr ~= sch_prev; % Trials switch
        tct = sch_curr == sch_prev; % Trials continuous

        tsw = logical([0 tsw]);
        tct = logical([0 tct]); % First trial is neither switch nor continuous - they're nothing!
        
        % Reduce switch trials to be indexed by samples
        tsw = tsw(md.sample_on_trials);
        tct = tct(md.sample_on_trials);
        
        % Display percent switch trials
        Ntsw = sum(tsw); Ntct = sum(tct);
        fprintf('Percent switch trials = %g \n', (Ntsw / (Ntsw + Ntct) *100) );
end


function X = convert_average(X,index)
% Converts X into an average in a specific frequency band
    X = X(index,:);
    X = mean(X);
    X = X(:);
end
