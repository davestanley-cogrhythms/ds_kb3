

function run_analysis_temp (sfc_mode, stage, file_range)

% % % Test Fries method on Kyle's simulated data. Then test Kyle's method
% and compare
% % % 
format compact
format long g

run_setdefaultfig
addpath(genpath('../funcs_supporting_local'));

if ~exist('sfc_mode','var'); sfc_mode = 2.2; end
    % sfc_mode = 1 -> fries
    % sfc_mode = 2 -> chronux
    % sfc_mode = 3 -> Lepage (not yet implemented)
    % sfc_mode = 4 -> mtcoherencysegc (segmented)
    % sfc_mode = 5 -> unit analysis
        % sfc_mode = 5.1 -> unit analysis, only this version reports the smoothed mean firing rates instad of firing rates per trial
    % sfc_mode = 6 -> raster_plotting
    % sfc_mode = *.2 -> use instad adjacent electrode! - not applicable for unit analysis (sfc_mode=5.*)
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
    file_range = [5];
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
    sfc_mode_group = floor(sfc_mode);
    sfc_mode_subgroup = round((sfc_mode - floor(sfc_mode))*10);
    exclude_clipping_trials=1;
    loadall_lfps = 1;      % Used for comparing the current spiking sequence against all LFP's
    
    
    if sfc_mode_subgroup == 2; do_adjacent = 1;
    else do_adjacent = 0; end
    
    if sfc_mode_group == 6     % Raster plotting mode
        curr_stage = 5; % For looking at the whole set of activity.
        plot_sfc_on = 1;
    end
    
    if sfc_mode == 5.1
        curr_stage = 5; % For looking at whole set of activity
    end
    
    apply_filter = 1;   % Apply sgolay filter to spiking data
    
    dt = get_dt;
    wind_s = round(0.2/dt);     % Window size for Fries method
    
    if do_adjacent && (sfc_mode ~= 5); adj_mode = '_adj'; else adj_mode = ''; end
    outpath = fullfile(getpath('path_buffer_curr'),['mode_' num2str(sfc_mode) adj_mode],['stage_' num2str(curr_stage)]);
    outname = fullfile(outpath,filename);
    if sfc_mode_group ~= 6; mkdir_dave(outpath); end
    if exist(outname,'file')
        fprintf('File %s exists...skipping\n',outname);
        %return
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
    
    N_categories = 4*2;
    ctgsetli = logical(zeros(length(get_good_samples(md.ctg1_trials,md)),N_categories+1));
    ctgsetli(:,1) = get_good_samples(md.ctg1_trials,md);
    ctgsetli(:,2) = get_good_samples(md.ctg2_trials,md);
    ctgsetli(:,3) = get_good_samples(md.ctg3_trials,md);
    ctgsetli(:,4) = get_good_samples(md.ctg4_trials,md);
    
    ctgsetli(:,5) = get_good_samples(md.ctg1_nr_trials,md);
    ctgsetli(:,6) = get_good_samples(md.ctg2_nr_trials,md);
    ctgsetli(:,7) = get_good_samples(md.ctg3_nr_trials,md);
    ctgsetli(:,8) = get_good_samples(md.ctg4_nr_trials,md);
    ctgsetli(:,9) = ones(size(ctgsetli,1),1);
    
    
    
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
        cat_range = 1:(N_categories+1);          % 9th category is all data!
        cat_range = 9;
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

            case 5
                sfc{i}.numspikes = [];
        end
        
        
        selected_units = 1:Nunits
        selected_units = 4;
        
        for j = selected_units
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
                switch sfc_mode_subgroup
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
                    sfc{i}.f = f;
                    sfc{i}.Nunits(j) = size(lfp_wind,2);
                    sfc{i}.spikerate_mu(j) = sum(sum(currspike_ctg,1)) / (size(currspike_ctg,1)*get_dt*sfc{i}.Ntraces(j));    % total # spikes / total time
                    %sfc{i}.spikes_per_trial{j} = sum(currspike_ctg,1)';
                    
                case 2
                    %tic
                    params.Fs = 1/dt;
                    params.tapers=get_tapers;

                    params.trialave = 1;
                    thresh = 5;
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

                    plot_sfc_all_electrodes = 0
                    if plot_sfc_all_electrodes
                        figure; plot(f,Cave);

                        figure
                        for jj = 1:size(currlfp2_ctg,3)
                            subplotsq(size(currlfp2_ctg,3),jj);
                            title([num2str(jj) '-' convert_unit_underscores(md.lfp_names(jj))])
                            [Cave3,phi2,S12ave2,S1ave2,S2ave2,f3]=coherencycpb(currlfp2_ctg(:,:,jj),currspike_ctg,params);
                            hold on; plot(f3,Cave3); ylim([0 0.25]);
                        end
                    end
                    
                    
                    % Plot individual trials and their statistics
                    params.trialave = 0;
                    [Cave_ind,phi,S12ave,S1ave,S2ave,f]=coherencycpb(currlfp_ctg,currspike_ctg,params);
                    freq_range = [ 15 30];
                    index = f >= freq_range(1) & f<=freq_range(2);
                    cave_stat = mean(Cave_ind(index,:),1);
                    figure;
%                     subplotsq(2,1); plott_matrix3D(f,Cave_ind,'showErrorbars',1); hold on; plot(f,mean(Cave_ind,2),'r','LineWidth',2); hold on; plot(f,Cave,'k','LineWidth',2)
                    subplotsq(2,1); h1 = plot(f,Cave_ind(:,1:min(end,5)),'m.'); hold on; h2 = plot(f,mean(Cave_ind,2),'r','LineWidth',2); hold on; h3 = plot(f,Cave,'k','LineWidth',2)
                    legend([h1(1),h2,h3],'SFC All Trials','Mean SFC','Estimator SFC');
                    subplotsq(2,2); hist(cave_stat,20)
                    title('SFC 14-30 Hz distribution')
                    xlabel('SFC'); ylabel('N');
                    
%                     for i = 1:size(currspike_ctg,2)
%                         [Cave_fries(:,i) f lfp_wind] = sfc_fries2D(currspike_ctg(:,i),currlfp_ctg(:,i),dt,wind_s);
%                         
%                         clf;
%                         subplot(10,1); plot(currlfp_ctg(:,i))
%                     end
                    

                    % Plot average sfc of "high" and "low" amplitude trials
                    sfc_high = cave_stat > 0.5;
                    sfc_low = ~sfc_high;
                    
                    params.trialave = 1;
                    [Cave_high,phi,S12ave_high,S1ave_high,S2ave_high,f]=coherencycpb(currlfp_ctg(:,sfc_high),currspike_ctg(:,sfc_high),params);
                    [Cave_low,phi,S12ave_low,S1ave_low,S2ave_low,f]=coherencycpb(currlfp_ctg(:,sfc_low),currspike_ctg(:,sfc_low),params);
                    figure; plot(f,Cave_high); hold on; plot(f,Cave_low,'r'); legend(['SFC upper N=' num2str(sum(sfc_high))],['SFC lower N=' num2str(sum(sfc_low))])
                    xlabel('Freq(hz)'); xlim([0 220]);
                    ylabel('SFC');
                    figure; plot(f,S1ave_high); hold on; plot(f,S1ave_low,'r'); legend('SFC upper','SFC lower')
%                     

                    THE CODE BELOW HERE IS GENERALL MESSED UP INDEX-WISE. NEED TO FIX

                    % Divide into ctgs into schemes that are continuous vs
                    % switch (a->b). 
                    clear sch
                    for iii = 1:Ntrials
                        [sch(iii)] = condition2morph(md.each_trial_condition(iii));
                    end
                    sch = (sch == 'A');
                    
                    sch_curr = sch(2:end);
                    sch_prev = sch(1:end-1);

                    tsw = sch_curr ~= sch_prev; % Trials switch
                    tct = sch_curr == sch_prev; % Trials cont
                    
                    tsw = logical([0 tsw]);
                    tct = logical([0 tct]);
                    
                    % PErcent of switch verus same trials
                    sum(tsw) / length(tsw) * 100
                    sum(tct) / length(tct) * 100
                    
                    
                    % Reduce dimensionality to only include good trials
                    switch_trial = tsw(md.sample_on_trials);
                    switch_trial = switch_trial(~curr_badtrials);
                    
                    same_trial = tct(md.sample_on_trials);
                    same_trial = same_trial(~curr_badtrials);
                    
                    sch = sch(md.sample_on_trials);
                    
                    keyboard
                    
%                     % Number fraction of scheme A's holding high vs low
%                     % trials
%                     This has an error in it - dimension mismatch between sch and sfc_high
%                     sum(sch(sfc_high)) / sum(sfc_high) * 100
%                     sum(sch(sfc_low)) / sum(sfc_low) * 100
                    
                    % Plot distribution of high sfc's
                    ctgsetli2 = ctgsetli(~curr_badtrials,:);
                    t = 1:size(ctgsetli2,1);
                    t = repmat(t(:),1,9);
                    figure; plott_matrix3D(t,ctgsetli2,'do_shift',2,'do_mean',0,'LineSpec',{'.'})
                    hold on; plott_matrix3D(t(sfc_high,:),ctgsetli2(sfc_high,:),'do_shift',2,'do_mean',0,'LineSpec',{'r.'})
                    
                    fract_high = [];
                    for ii = 1:9
                        fract_high(ii) = sum(ctgsetli2(sfc_high,ii)) / sum(ctgsetli2(:,ii))*100
                    end
                    
                    
%                     clear numspikes
%                     for  ii = 1:9
%                         numspikes_temp = (get_numspikes_fromfile(curr_unit,curr_stage,md));
%                         %                         % Normalizing
%                         %                         numspikes = (numspikes - min(numspikes)) / (max(numspikes) - min(numspikes));
%                         numspikes{ii} = (numspikes_temp(ctgsetli(:,ii)));
%                         muspikes(ii) = mean(numspikes{ii});
%                     end
%                     
%                     figure; plot(fract_high, muspikes,'b.')
%                     for ii = 1:9
%                         hold on; text(fract_high(i), muspikes(i),['' num2str(i)])
%                     end
%                     xlabel('Ctg Fract high beta');
%                     ylabel('Num spikes');
                    
%                     % If the chosen set of trials contains trial 1, we need to drop it, since there is no previous trial. We then restore the vector to the correct length below by adding zero to both, so the 1st trial is removed from the average
%                     % Drop first trial. This code is BROKEN
%                     
%                     ctgsetli2 = ctgsetli(~curr_badtrials,:);
%                     ctg_sets = find(ctgsetli2(:,i));
% 
%                     drop_first = 0;
%                     if ctg_sets(1) == 1
%                         ctg_sets = ctg_sets(2:end);
%                         drop_first = 1;
%                     end
%                     sch_curr = sch(ctg_sets);
%                     sch_prev = sch(ctg_sets-1);
%                     
%                     same_trial = sch_curr == sch_prev;
%                     switch_trial = sch_curr ~= sch_prev;
%                     if drop_first
%                         same_trial = logical([0 same_trial]);
%                         switch_trial = logical([0 switch_trial]);
%                     end
                    
                    params.trialave = 1;
                    [Cave_same,phi,S12ave_high,S1ave_high,S2ave_high,f]=coherencycpb(currlfp_ctg(:,same_trial),currspike_ctg(:,same_trial),params);
                    [Cave_switch,phi,S12ave_low,S1ave_low,S2ave_low,f]=coherencycpb(currlfp_ctg(:,switch_trial),currspike_ctg(:,switch_trial),params);
                    
                    figure; plot(f,Cave_same); hold on; plot(f,Cave_switch,'r'); legend(['same N=' num2str(sum(same_trial))],['switch N=' num2str(sum(switch_trial))])

                    
                    
                    % Compare each trial's response to sfc value
                    etr = md.each_trial_response;
                    etr = etr(md.sample_on_trials);
                    etr = etr(~curr_badtrials);
                    
                    
                    
% 
%                     schA = any(ctgsetli(:,1:2),2);
%                     schB = any(ctgsetli(:,3:4),2);
%                     figure; imagesc([schA schB])
                    
                    
                    keyboard
                    
                    sfc{i}.Cave = [sfc{i}.Cave Cave(:)];
                    sfc{i}.S12ave = [sfc{i}.S12ave S12ave(:)];
                    sfc{i}.S1ave = [sfc{i}.S1ave S1ave(:)];
                    sfc{i}.S2ave = [sfc{i}.S2ave S2ave(:)];
                    sfc{i}.f = f;
                    
                    
                    
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
        %save(outname,'sfc');
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