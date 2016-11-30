

function run_stats_specialization (sfc_mode,stage_range, file_range)

    % % % Test Fries method on Kyle's simulated data. Then test Kyle's method
    % and compare
    % % % 
    format compact
    format long g

    run_setdefaultfig
    if ~isdeployed
        addpath(genpath('./funcs_supporting_local'));
    end

    if ~exist('sfc_mode','var'); sfc_mode = 5; end
        % sfc_mode = 5 -> Calculates raw firing rates in each ctg. Some of
        % these are morph specific, so it is possible to calculte the
        % categorizaiton index as well

    if ~exist('stage','var'); stage_range = [2,3]; end
        % stage = 1 - pre sample on
        % stage = 2 - sample stage
        % stage = 3 - delay stage
        % stage = 4 - sample + delay stage
        % stage = 5 - all
        % stage = 6 - Test1 - match or non-match
        % stage = 7 - TDelay - test delay
        % stage = 8 - Test2 - match (if prior was non-match)

   
    filepaths = fullfile(getpath('path_metadata'));
    file_list = get_filelist(filepaths);
    Nfiles = length(file_list);


    if  ~exist('file_range','var')
       %file_range = 1:8;
       file_range = [1:Nfiles];
    end

    for curr_stage = stage_range
        stats.(stagename(curr_stage)) = calc_specialization_files(filepaths,file_list,file_range,sfc_mode,curr_stage);
    end
    
    path_out = getpath('path_buffer_specialization_stats');
    
    %dont forget to update getstages ir
    
    mkdir_dave(path_out)
    save(fullfile(path_out,'specialization.mat'),'stats');
    
end

function stats_stage = calc_specialization_files(filepaths,file_list,file_range,sfc_mode,curr_stage)

    parfor i = file_range
        fprintf ('Processing file %d of %d, name %s \n',i,file_range(end),file_list{i});
        stats(i) = calc_specialization(filepaths,file_list{i},sfc_mode,curr_stage);
    end
    
    %save temp.mat
    
    mu_arr = [];
    std_arr = [];
    for i = file_range
        temp_mu=cat(2,stats(i).data{:});
        temp_mu = arrayfunu(@(s) mean(s.spikerates),temp_mu);
        temp_mu = vertcat(temp_mu{:});
        temp_mu = permute(temp_mu,[3,2,1]);
        mu_arr = cat(2,mu_arr,temp_mu);
        
        temp_std=cat(2,stats(i).data{:});
        temp_std = arrayfunu(@(s) mean(s.spikerates),temp_std);
        temp_std = vertcat(temp_std{:});
        temp_std = permute(temp_std,[3,2,1]);
        std_arr = cat(2,std_arr,temp_std);
    end
    
    stats_stage.stats = stats;
    stats_stage.mu_arr = mu_arr;
    stats_stage.std_arr = std_arr;
    

end


function stats = calc_specialization(filepaths,filename_orig,sfc_mode,curr_stage)
    

    % Generate outpath
    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    sfc_mode_group = floor(sfc_mode); 

    dt = get_dt;
    
    % Function specific options
    plot_on = 0;
    
    md = load (fullfile(filepaths,filename_orig));     % Load metadata
    
    % Build metadata structure containing parameters for this function
    Nsamples = length(md.sample_on_times); %Number of "good" samples - when trial reached sample_on stage
    Nunits = length(md.unit_names);    % Number of units tracked
    Nelects = length(md.lfp_names);      % Number of electrodes
    
    % Setup ctgsetli
    [ctgsetli,ctgsetnames] = get_ctgsetli_sp(md);
    N_categories = size(ctgsetli,2);    % Update N_categories
    
    % Load data appropriate to current stage
    [currspike_all] = get_spike_ts_all(Nunits,curr_stage,md);
    
    ticID = tic;
    for i = 1:N_categories
        

        data{i}.spikerates = [];
        
        for j = 1:Nunits

            
            [i j];
            curr_unit = j;
            unit_names2= strrep(md.unit_names,'_',' ');
            
            switch mode_subgroups(1)    % Check if mode is to do unit analysis trace plotting
                case 0
                    spikerates = get_numspikes_fromfile(curr_unit,curr_stage,md);
%                         % Normalizing
%                         numspikes = (numspikes - min(numspikes)) / (max(numspikes) - min(numspikes));
                    spikerates = spikerates(ctgsetli(:,i));
                case 1
                    [currspike1] = get_spike_ts(curr_unit,curr_stage,md);
                    currspike1 = currspike1(:,ctgsetli(:,i));
                    spikerates = mean(currspike1,2); spikerates = spikerates / dt;

                    apply_filter = 1;   % Apply sgolay filter to spiking data
                    if apply_filter
                        spikerates = sgolayfilt(spikerates,3,round(0.151/dt));   % Apply 3rd-order filter, 151 ms 
                    end
            end
                        
            switch sfc_mode_group
                case 5
                    data{i}.spikerates = [data{i}.spikerates spikerates(:)];
            end
            
        end

    end
    elapsed_time = toc(ticID);
    
    % Build metadata structure containing parameters for this function
    m.sfc_mode = sfc_mode;
    m.curr_stage = curr_stage;
    m.Nsamples = Nsamples; %Number of "good" samples - when trial reached sample_on stage
    m.Nunits = Nunits;    % Number of units tracked
    m.Nelects = Nelects;      % Number of electrodes
    m.ctgsetli=ctgsetli;
    m.ctgsetnames=ctgsetnames;
    m.elapsed_time = elapsed_time;
    m.filename = filename_orig;
    m.unit_names = md.unit_names;
    
    stats.m = m;
    stats.data = data;

end





function [ctgsetli,ctgsetnames] = get_ctgsetli_sp(md)

    test_on = 0;
    plot_on = 0;
    
    Ntrials = length(md.start_times);
    for iii = 1:Ntrials
        [sch(iii), morphA(iii), morphB(iii)] = condition2morph(md.each_trial_condition(iii)); % Slow, but not limiting.
    end
    
    % Test nomenclature for A1,A2,B1,B2
    if test_on
        A1 = get_good_samples(md.ctg1_trials,md);
        A1_test = get_good_samples(find(morphA(:) > 50 & morphB(:) ~= 50 & sch(:)=='A') ,md);
        
        A2 = get_good_samples(md.ctg2_trials,md);
        A2_test = get_good_samples(find(morphA(:) < 50 & morphB(:) ~= 50 & sch(:)=='A') ,md);
        
        B1 = get_good_samples(md.ctg3_trials,md);
        B1_test = get_good_samples(find(morphB(:) > 50 & morphA(:) ~= 50 & sch(:)=='B') ,md);
        
        B2 = get_good_samples(md.ctg4_trials,md);
        B2_test = get_good_samples(find(morphB(:) < 50 & morphA(:) ~= 50 & sch(:)=='B') ,md);
        
        sum(A1 ~= A1_test)
        sum(A2 ~= A2_test)
        sum(B1 ~= B1_test)
        sum(B2 ~= B2_test)
    end
    
    clear ctgsetli ctgsetnames
    % Scheme A specializations
    i=0;
    i=i+1; ctgsetnames{i} = 'Sp_SchA_A1B1'; ctgsetli_tr(:,i) = sch(:)=='A' & morphA(:) > 50 & morphB(:) > 50;
    i=i+1; ctgsetnames{i} = 'Sp_SchA_A1B2'; ctgsetli_tr(:,i) = sch(:)=='A' & morphA(:) > 50 & morphB(:) < 50;
    i=i+1; ctgsetnames{i} = 'Sp_SchA_A2B1'; ctgsetli_tr(:,i) = sch(:)=='A' & morphA(:) < 50 & morphB(:) > 50;
    i=i+1; ctgsetnames{i} = 'Sp_SchA_A2B2'; ctgsetli_tr(:,i) = sch(:)=='A' & morphA(:) < 50 & morphB(:) < 50;
    
    % Scheme B specializations
    i=i+1; ctgsetnames{i} = 'Sp_SchB_A1B1'; ctgsetli_tr(:,i) = sch(:)=='B' & morphA(:) > 50 & morphB(:) > 50;
    i=i+1; ctgsetnames{i} = 'Sp_SchB_A1B2'; ctgsetli_tr(:,i) = sch(:)=='B' & morphA(:) > 50 & morphB(:) < 50;
    i=i+1; ctgsetnames{i} = 'Sp_SchB_A2B1'; ctgsetli_tr(:,i) = sch(:)=='B' & morphA(:) < 50 & morphB(:) > 50;
    i=i+1; ctgsetnames{i} = 'Sp_SchB_A2B2'; ctgsetli_tr(:,i) = sch(:)=='B' & morphA(:) < 50 & morphB(:) < 50;
    
    % Scheme A rEL categorizations
    i=i+1; ctgsetnames{i} = 'Ca_SchA_MA100'; ctgsetli_tr(:,i) = sch(:)=='A' & morphA(:) == 100 & morphB(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchA_MA080'; ctgsetli_tr(:,i) = sch(:)=='A' & morphA(:) ==  80 & morphB(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchA_MA060'; ctgsetli_tr(:,i) = sch(:)=='A' & morphA(:) ==  60 & morphB(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchA_MA040'; ctgsetli_tr(:,i) = sch(:)=='A' & morphA(:) ==  40 & morphB(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchA_MA020'; ctgsetli_tr(:,i) = sch(:)=='A' & morphA(:) ==  20 & morphB(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchA_MA000'; ctgsetli_tr(:,i) = sch(:)=='A' & morphA(:) ==   0 & morphB(:) ~= 50;
    
    % Scheme A NR categorizations
    i=i+1; ctgsetnames{i} = 'Ca_SchA_MB100'; ctgsetli_tr(:,i) = sch(:)=='A' & morphB(:) == 100 & morphA(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchA_MB080'; ctgsetli_tr(:,i) = sch(:)=='A' & morphB(:) ==  80 & morphA(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchA_MB060'; ctgsetli_tr(:,i) = sch(:)=='A' & morphB(:) ==  60 & morphA(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchA_MB040'; ctgsetli_tr(:,i) = sch(:)=='A' & morphB(:) ==  40 & morphA(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchA_MB020'; ctgsetli_tr(:,i) = sch(:)=='A' & morphB(:) ==  20 & morphA(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchA_MB000'; ctgsetli_tr(:,i) = sch(:)=='A' & morphB(:) ==   0 & morphA(:) ~= 50;
    
    % Scheme B REL categorizations
    i=i+1; ctgsetnames{i} = 'Ca_SchB_MB100'; ctgsetli_tr(:,i) = sch(:)=='B' & morphB(:) == 100 & morphA(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchB_MB080'; ctgsetli_tr(:,i) = sch(:)=='B' & morphB(:) ==  80 & morphA(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchB_MB060'; ctgsetli_tr(:,i) = sch(:)=='B' & morphB(:) ==  60 & morphA(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchB_MB040'; ctgsetli_tr(:,i) = sch(:)=='B' & morphB(:) ==  40 & morphA(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchB_MB020'; ctgsetli_tr(:,i) = sch(:)=='B' & morphB(:) ==  20 & morphA(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchB_MB000'; ctgsetli_tr(:,i) = sch(:)=='B' & morphB(:) ==   0 & morphA(:) ~= 50;
    
    % Scheme B NR categorizations
    i=i+1; ctgsetnames{i} = 'Ca_SchB_MA100'; ctgsetli_tr(:,i) = sch(:)=='B' & morphA(:) == 100 & morphB(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchB_MA080'; ctgsetli_tr(:,i) = sch(:)=='B' & morphA(:) ==  80 & morphB(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchB_MA060'; ctgsetli_tr(:,i) = sch(:)=='B' & morphA(:) ==  60 & morphB(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchB_MA040'; ctgsetli_tr(:,i) = sch(:)=='B' & morphA(:) ==  40 & morphB(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchB_MA020'; ctgsetli_tr(:,i) = sch(:)=='B' & morphA(:) ==  20 & morphB(:) ~= 50;
    i=i+1; ctgsetnames{i} = 'Ca_SchB_MA000'; ctgsetli_tr(:,i) = sch(:)=='B' & morphA(:) ==   0 & morphB(:) ~= 50;
    
    % All
    i=i+1; ctgsetnames{i} = 'Al_Good_00000'; ctgsetli_tr(:,i) = true(Ntrials,1);
    
    % Calculate correct trials
    Nctgs = size(ctgsetli_tr,2);
    corrtrcd = 0;
    correct_trials = md.each_trial_response(:) == corrtrcd;
    
    % Include only correct trials
    if plot_on; figure; bar(sum(ctgsetli_tr)); end
    ctgsetli_tr = ctgsetli_tr & repmat(correct_trials,1,Nctgs);
    if plot_on; hold on; bar(sum(ctgsetli_tr),'r'); xlabel('Ctgset'); ylabel('Number of trials'); legend('All trials','All good');end
    
    % Convert trials to samples
    Nsamples = sum(md.sample_on_trials);
    ctgsetli = false(Nsamples, Nctgs);
    for i = 1:Nctgs
            chosen_trials = ctgsetli_tr(:,i);
            samplesli = (trials_to_samples( chosen_trials, md.sample_on_trials));
            ctgsetli(:,i) = [samplesli(:)];
    end

end



function currspike_all = get_spike_ts_all(Nunits,curr_stage,md)


    j = 1;
    [currspike] = get_spike_ts(j,curr_stage,md);
    sz = size(currspike);
    currspike_all = zeros([sz, Nunits]);
    currspike_all(:,:,j) = currspike;
    
    for j = 1:Nunits
        currspike_all(:,:,j) = get_spike_ts(j,curr_stage,md);
    end
    
end

function y = do_baseline_subtract(x,baseline_subtract)
    if baseline_subtract
        sz = size(x);
        y = x - repmat(mean(x,2),1,sz(2));
    else
        y = x;
    end
end


function outname = stagename(curr_stage)


    outname = num2str(curr_stage);
    outname = strrep(outname,'-','m');
    outname = ['stage' outname];

end