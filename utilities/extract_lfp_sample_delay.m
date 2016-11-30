
function extract_lfp_sample_delay

    outpath = fullfile(getpath('path_lfp_sample_delay'));
    
    mkdir_dave (outpath);
    
    file_list = dir(fullfile(getpath('path_lfp'),'*.mat'));
    
    parfor i = 1:length(file_list)
        fprintf ('Processing file %d of %d, name %s \n',i,length(file_list),file_list(i).name);
        extract_sample_delay(file_list(i).name,outpath);
    end
    
end


function extract_sample_delay(filename, outpath)


    % Note: indices refers to indices of lfp matrix
    % Trials refers to indices of the ctgX_trials arrays
    
    plot_on = 0;
    
    if ~exist('filename','var'); filename = 'L091906.mat'; end
    metadata_lfp = [];
    
    
    if exist(fullfile(outpath,filename),'file')
        fprintf('File already exists. Skipping.\n');
        return;
    end

    load (fullfile(getpath('path_metadata'),filename));     % Load metadata
    load (fullfile(getpath('path_lfp'),filename));     % Load full file
    lfpmat = double(lfpmat_int16); clear lfpmat_int16;
    
    
    format long g
    t = (0:size(lfpmat,1)-1) * dt;
    stagesir_all = get_stagesir(5);
    Nsd = stagesir_all(2)+1;  % sample duration - Number of indices in 1700 msec - allow 100ms for transmission delay to PFC
    Nsd_pre = abs(stagesir_all(1));

    Ntrials = length(start_times);  % Total number of trials
    Nsamples = length(sample_on_times); %Number of "good" samples - when trial reached sample_on stage
    Nunits = length(unit_names);    % Number of units tracked
    Nelects = length(lfp_names);      % Number of electrodes

    
    metadata_lfp.good_samples = (sample_on_ind + Nsd - 1 < length(lfpmat) & sample_on_ind - Nsd_pre > 0);
    metadata_lfp.good_trials = ss_ind(:,2) < length(lfpmat);
    
    if sum(metadata_lfp.good_samples) < length(sample_on_times)
        
        fprintf('**Warning!! More samples than LFP! Trimming...**\n');
        %pause;
        metadata_lfp.samples_dropped = 1;
        sample_on_ind2 = sample_on_ind( metadata_lfp.good_samples );
        Nsamples = length(sample_on_ind2);
        %return;
        
    else
        metadata_lfp.samples_dropped = 0;
        sample_on_ind2 = sample_on_ind;
    end
    
    
    % Pull out all sample_ons.
    sample_and_delay_indices = repmat(sample_on_ind2(:)',(Nsd+Nsd_pre),1) + repmat(( -Nsd_pre:(Nsd-1)  )',1,Nsamples);
    sample_and_delay_indices = sample_and_delay_indices(:);
    lfp_samples = lfpmat(sample_and_delay_indices,:);
    lfp_samples = permute(lfp_samples,[1 3 2]);
    lfp_samples = reshape(lfp_samples,Nsd+Nsd_pre,Nsamples,Nelects);
    %figure; plot(t,lfpmat(:,1)); hold on;
    %plot(t(sample_and_delay_indices(:)),lfpmat(sample_and_delay_indices(:),1),'g')
    
    lfp_sample_indices = [sample_on_ind2(:)-Nsd_pre sample_on_ind2(:) sample_on_ind2(:)+Nsd-1]; % Start_index Sample on End_index
    
    % Format of lfp_sample is lfp_sample(sample data, sample number, electrode number)
    lfp_sample = int16(lfp_samples);
    save(fullfile(outpath,filename),'lfp_sample','metadata_lfp','lfp_sample_indices');
    
    if plot_on
        shortening_factor = 10;
        enum = 1;
        figure;
        plot(lfpmat(1:round(length(lfpmat)/shortening_factor),enum));
        
        
        N_samples_to_plot= round(Nsamples / shortening_factor);
        hold on;
        plotting_indices = [];
        for j = 1:N_samples_to_plot
            plotting_indices = [plotting_indices (lfp_sample_indices(j,1):lfp_sample_indices(j,3))'];
        end
        plot( plotting_indices,lfp_sample(:,1:N_samples_to_plot,enum));
        
    end
 
    clear lfpmat
end


