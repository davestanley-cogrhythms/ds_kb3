

function calc_coherence(filename, sfc_mode)

    if ~exist('sfc_mode','var'); sfc_mode = 4; end
                % sfc_mode = 1 -> fries
                % sfc_mode = 2 -> chronux
                % sfc_mode = 3 -> Lepage (not yet implemented)
                % sfc_mode = 4 -> mtcoherencysegc (segmented)

    % Note: indices refers to indices of lfp matrix
    % Trials refers to indices of the ctgX_trials arrays
    
    if ~exist('filename','var'); filename = 'L091906.mat'; end

    load (fullfile(getpath('path_metadata'),filename));     % Load metadata
    %load (fullfile(getpath('path_lfp'),filename));     % Load full file
    %lfpmat = double(lfpmat_int16); clear lfpmat_int16;
    
    load (fullfile(getpath('path_lfp_sample_delay'),filename));
    lfp_sample = double(lfp_sample);
    
    
    
    plot_on = 0;
    plot_debug = 0;
    run_tests = 0;
    
    
    format long g
    %t = (0:size(lfpmat,1)-1) * dt;
    
    Ntrials = length(start_times);  % Total number of trials
    Nsamples = length(sample_on_times); %Number of "good" samples - when trial reached sample_on stage
    Nunits = length(unit_names);    % Number of units tracked
    Nelects = length(lfp_names);      % Number of electrodes

    unum = 1;   % Unit number to plot
    enum = 1;   % Electrode number to plot

    % % Do some basic sanity checks -- some things don't make sense here.
    % Need to come back and revisit.
    if run_tests
        
        Nat = length(start_times);
        
        % First, all "correct" trials should be a subset of "sample_on_trials"
        out = is_subset(sample_on_trials(:),each_trial_response == 0)
        figure; plot(each_trial_response == 0,'r.')
        hold on; plot(sample_on_trials,'.')
            % All "correct" are contained within the "sample_on_trials"
        
            
        % Second, all correct trials should be a subset of ctgx_trials?
        all_trials = ones(1,Nat);
        all_ctgX = (build_logical(ctg1_trials,Nat) == all_trials) | (build_logical(ctg2_trials,Nat) == all_trials) ...
            | (build_logical(ctg3_trials, Nat) == all_trials) | (build_logical(ctg4_trials, Nat) == all_trials);
        
        out = is_subset(all_ctgX(:), each_trial_response == 0)
            %No also?! Are there outliers due to the ones that lie on the
            %border?
            
        % Third, all ctgx_trials should be a subset of "sample_on_trials"
        Nat = length(start_times);
        all_trials = ones(1,Nat);
        all_ctgX = (build_logical(ctg1_trials,Nat) == all_trials) | (build_logical(ctg2_trials,Nat) == all_trials) ...
            | (build_logical(ctg3_trials, Nat) == all_trials) | (build_logical(ctg4_trials, Nat) == all_trials);
        
        out = is_subset(sample_on_trials,all_ctgX) 
            % Example - point 230
        each_trial_response(230)
        find(ctg1_trials == 230)
        find(ctg2_trials == 230)
        find(ctg3_trials == 230)
        find(ctg4_trials == 230)
        find(ctg1_nr_trials == 230)
        find(ctg2_nr_trials == 230)
        find(ctg3_nr_trials == 230)
        find(ctg4_nr_trials == 230)

        figure; plot(sample_on_trials,'.')
        figure; plot(all_ctgX,'g.')
            % Some of the ctgX trials are outside of the set of
            % sample_on_trials. I guess this means that these trials
            % failed.
         
    end
    
    
    
    % Don't need to do this since I already have this data in saved files
%     % Pull out all sample_ons.
%     Nsd = round(1.6 / dt);  % sample duration - Number of indices in 1600 msec
%     sample_and_delay_indices = repmat(sample_on_ind(:)',Nsd,1) + repmat((1:Nsd)',1,Nsamples);
%     sample_and_delay_indices = sample_and_delay_indices(:);
%     lfp_sample = lfpmat(sample_and_delay_indices,:);
%     lfp_sample = permute(lfp_sample,[1 3 2]);
%     lfp_sample = reshape(lfp_sample,Nsd,Nsamples,Nelects);
%     %figure; plot(t,lfpmat(:,1)); hold on;
%     %plot(t(sample_and_delay_indices(:)),lfpmat(sample_and_delay_indices(:),1),'g')
%     %clear lfpmat

    % Pull out all trials (we won't be needing these)
%     Ntd_all = diff(ss_ind,[],2);
%     Ntd = round(max([median(Ntd_all) sfc_mode(Ntd_all) mean(Ntd_all)]));  % trial duration in indices
%     sample_and_delay_indices = repmat(ss_ind(:,1)',Ntd,1) + repmat((1:Ntd)',1,Ntrials);
%     sample_and_delay_indices = sample_and_delay_indices(:);
%     lfp_trials = lfpmat(sample_and_delay_indices,:);
%     lfp_trials = permute(lfp_trials,[1 3 2]);
%     lfp_trials = reshape(lfp_trials,Ntd,Ntrials,Nelects);
%     figure; plot(t,lfpmat(:,1)); hold on; plot(t(sample_and_delay_indices(:)),lfpmat(sample_and_delay_indices(:),1),'g')


%     figure;
%     plot(t,lfp_matrix(:,enum));
%     hold on;
%     plot(t(unit_ind{unum}),lfp_matrix(unit_ind{unum},enum),'r.');

    % Test run to get dimensionality
    i=1;j=1;
    Nsd = size(lfp_sample,1);
    curr_electrode = unit_LFP_pairs(2,i);
    currlfp = (lfp_sample(:,i,curr_electrode));
    spike_trace = get_spikes_n_range(unit_ind{i},[sample_on_ind(j) (sample_on_ind(j) + Nsd - 1)]);
    trange = (0:Nsd-1)*dt;


    switch sfc_mode
        case 1
            [C f] = sfc_fries([spike_trace(:) currlfp(:)],dt,round(0.5/dt));
            C_arr = zeros(length(C),Nsamples,Nunits);
        case 2
            params.Fs = 1/dt;
            params.tapers=[3 5];
            [C,phi,S12,S1,S2,f]=coherencycpb(currlfp(:),spike_trace(:),params);
        case 4
            params.Fs = 1/dt;
            params.tapers=[3 5];
            segave=0;
            win = [0.2];
            [C,phi,S12,S1,S2,f,zerosp]=coherencysegcpb(currlfp(:),spike_trace(:),win,params,segave);

    end

    
    i=1;
    for i = 1:Nunits
       
        for j = 1:Nsamples
            curr_electrode = unit_LFP_pairs(2,i);
            currlfp = (lfp_sample(:,j,curr_electrode));
            spike_trace = get_spikes_n_range(unit_ind{i},[sample_on_ind(j) (sample_on_ind(j) + Nsd - 1)]);
            trange = (0:Nsd-1)*dt;
            j;
            if plot_on figure; plot(trange,zscore([spike_trace(:) currlfp(:)],[])); end
            switch sfc_mode
                case 1
                    [C ftemp] = sfc_fries([spike_trace(:) currlfp(:)],dt,round(0.5/dt));
                    C_arr(:,j,i) = C(:);
                    if sum(isnan(ftemp)) == 0; f = ftemp; end   % ftemp may be all NaNs. Therefore, only save it if it's not.
                case 2
                    params.Fs = 1/dt;
                    params.tapers=[3 5];
                    [C,phi,S12,S1,S2,f]=coherencycpb(currlfp(:),spike_trace(:),params);
                    C_arr(:,j,i) = C(:);
                    phi_arr(:,j,i) = phi(:);
                    S12_arr(:,j,i) = S12(:);
                    S1_arr(:,j,i) = S1(:);
                    S2_arr(:,j,i) = S2(:);
                    
                    
            end
            if plot_on

            end

        end
    end
    

    
    if plot_on 
        alpha = [8 12];
        beta = [13 30];
        gamma = [30 80];

        aind = (f >= alpha(1)) & (f < alpha(2));
        bind = (f >= beta(1)) & (f < beta(2));
        gind = (f >= gamma(1)) & (f < gamma(2));

        acoh = mean(C_arr(aind,:,:));
        bcoh = mean(C_arr(bind,:,:));
        gcoh = mean(C_arr(gind,:,:));
        
        % % This is borkenheid!
        unit_to_plot = 1;
        good_ind = sum(isnan(C_arr)) == 0;
        C_arr_clean = C_arr(:,good_ind(:,:,unit_to_plot),unit_to_plot);
        figure; plot(f, mean(C_arr_clean,2));
%         figure; plot(f, mean(S1,2));
%         figure; plot(f, mean(S2,2));
        
        acoh = mean(C_arr_clean(aind,:,:));
        bcoh = mean(C_arr_clean(bind,:,:));
        gcoh = mean(C_arr_clean(gind,:,:));
        
        % % Plot individual trials
        unit_to_plot = 1;
        good_ind = sum(isnan(C_arr(:,:,unit_to_plot)),1) == 0;
        chosen_ctg_trials = build_logical(ctg1_trials,Ntrials); tr1_sind = (trials_to_samples( chosen_ctg_trials(:) & each_trial_response(:) == 0, sample_on_trials) & good_ind);
        chosen_ctg_trials = build_logical(ctg2_trials,Ntrials); tr2_sind = (trials_to_samples( chosen_ctg_trials(:) & each_trial_response(:) == 0, sample_on_trials) & good_ind);
        chosen_ctg_trials = build_logical(ctg3_trials,Ntrials); tr3_sind = (trials_to_samples( chosen_ctg_trials(:) & each_trial_response(:) == 0, sample_on_trials) & good_ind);
        chosen_ctg_trials = build_logical(ctg4_trials,Ntrials); tr4_sind = (trials_to_samples( chosen_ctg_trials(:) & each_trial_response(:) == 0, sample_on_trials) & good_ind);
        
        
        mat = C_arr(:,:,unit_to_plot); % Matrix to plot
        good_ind = sum(isnan(mat(:,:)),1) == 0;
        
        figure;
        hold on; sinds = tr1_sind; errorbar_dave(f,mean(mat(:,sinds),2),std(mat(:,sinds),[],2),'b');
        hold on; sinds = tr2_sind; errorbar_dave(f,mean(mat(:,sinds),2),std(mat(:,sinds),[],2),'g');
        hold on; sinds = tr3_sind; errorbar_dave(f,mean(mat(:,sinds),2),std(mat(:,sinds),[],2),'r');
        hold on; sinds = tr4_sind; errorbar_dave(f,mean(mat(:,sinds),2),std(mat(:,sinds),[],2),'m');
        legend('SchA 1','SchA 2','SchB 1','SchB 2')
        xlabel('freq (Hz)');ylabel('coherence');
        title(filename);
        
        figure;
        hold on; sinds = tr1_sind; plot(f,mean(mat(:,sinds),2),'b');
        hold on; sinds = tr2_sind; plot(f,mean(mat(:,sinds),2),'g');
        hold on; sinds = tr3_sind; plot(f,mean(mat(:,sinds),2),'r');
        hold on; sinds = tr4_sind; plot(f,mean(mat(:,sinds),2),'m');
        legend('SchA 1','SchA 2','SchB 1','SchB 2')
        xlabel('freq (Hz)');ylabel('coherence');
        title(filename);
        
        
        % % Calc and plot the SFC for the whole dataset. (requires super
        % computer!)
        
        %temp = load (fullfile(getpath('path_roydata_orig'),filename));     % Load full file

        load (fullfile(getpath('path_lfp'),filename));     % Load full file
        lfpmat = double(lfpmat_int16); clear lfpmat_int16;
        
        i=1;j=1;
        Nsd = length(lfpmat);
        curr_electrode = unit_LFP_pairs(2,i);
        currlfp = lfpmat(:,curr_electrode);
        spike_trace = get_spikes_n_range(unit_ind{i},[1 length(lfpmat)]);
        trange = (0:Nsd-1)*dt;


        switch sfc_mode
            case 1
                [C f] = sfc_fries([spike_trace(:) currlfp(:)],dt,round(0.5/dt));

            case 2
                params.Fs = 1/dt;
                params.tapers=[3 5];
                [C_sm,phi,S12,S1,S2,f]=coherencycpb(currlfp(:),spike_trace(:),params);

           

                df = mode(diff(f));
                Nsmooth = round(0.5 / df);
                C_sm = wkeep(conv(C,1/Nsmooth*ones(1,Nsmooth)),length(C),'c');

        end
        figure; plot(f,C_sm);
        title(filename);
        xlabel('freq (Hz)');ylabel('coherence');
 
    end
    
    [pathstr, name, ext] = fileparts(filename);
    %mkdir_dave(getpath('path_SFC'));
    
    switch sfc_mode
        case 1
            save(fullfile(getpath('path_SFC'),[name '_sfc_mode_' num2str(sfc_mode) '.mat']),'C_arr','f');
        case 2
            save(fullfile(getpath('path_SFC'),[name '_sfc_mode_' num2str(sfc_mode) '.mat']),'C_arr','S1_arr','f');
            
    end
%     figure;
%     ind = 42;
%     tli = ss_ind(ind,1):ss_ind(ind,2);
%     plot(t(tli),lfp_matrix(tli,enum),'k');
%     tui = unit_ind{unum}( unit_ind{unum} >= ss_ind(ind,1) & unit_ind{unum} <= ss_ind(ind,2)  );
% 
%     hold on; plot(t(tui),lfp_matrix(tui,enum),'m.')

end





