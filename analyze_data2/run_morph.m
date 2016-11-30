

function run_morph(stage,morph_mode)

    format compact
    format long g

    run_setdefaultfig
    addpath(genpath('./funcs_supporting_local'));

    if ~exist('morph_mode','var'); morph_mode = 2; end
        % So far there is only one morph mode

    if ~exist('stage','var'); stage = 3; end
        % stage = 1 - pre sample on
        % stage = 2 - sample stage
        % stage = 3 - delay stage
        % stage = 4 - sample + delay stage
        % stage = 5 - all


    % path_metadata = getpath('path_metadata');
    

    bad_trials = prep_badtrials(stage);


    file_list = get_filelist(getpath('path_metadata'));
    Nfiles = length(file_list);
    
    outpath = fullfile(getpath('path_buffer_curr'),'morphs');
    filename = ['mode_' num2str(morph_mode) '_stage_' num2str(stage)];
    mkdir_dave(outpath)
    outname = fullfile(outpath,filename);

    if  ~exist('file_range','var')
        file_range = 1:Nfiles;
        %file_range = [1:5];
    end

    parfor i = file_range
        fprintf ('Processing file %d of %d, name %s \n',i,length(file_list),file_list{i});
        [sens{i},unit_names{i}] = calc_stats(file_list{i},morph_mode,stage,bad_trials,i);
    end
    fnames = file_list(file_range);
    funames = fname_uname_merge (fnames,unit_names);  % File and lfp names grouped by file
    
    save(outname,'sens','file_list','funames');
    
    sens_merg = vertcat(sens{:});


end


function [sens,unit_names] = calc_stats(filename,morph_mode,stage,bad_trials,i)

    debug_on = 0;
    %plot_on = 1;

    % Load metadata
    md = load(fullfile(getpath('path_metadata'),filename));

    
    % Get dimensions & setup parameters
    Nunits = length(md.unit_names);
    Ntrials = length(md.sample_on_trials);
    Nsamples = sum(md.sample_on_trials);
    corrtrcd = 0; % Code for correct response
    
    % Import data from structure
    sample_on_trials = md.sample_on_trials;
    each_trial_condition = md.each_trial_condition;
    each_trial_response = md.each_trial_response;
    correct_trials = each_trial_response(:) == corrtrcd;
    
    % Samples subset
    correct_samples = correct_trials(sample_on_trials);
    each_sample_condition = each_trial_condition(sample_on_trials);
    
    % Correct samples subset
    each_cs_condition = each_sample_condition(correct_samples);     % Conditions for only Correct Samples (CS)
    Ncs = length(each_cs_condition);
    
    morphsA = zeros(1,Ncs);
    morphsB = zeros(1,Ncs);
    sch = char(zeros(1,Ncs));
    for i = 1:Ncs
        [sch(i) morphA(i) morphB(i)] = condition2morph(each_cs_condition(i));
    end
    
    
    
    
    % Total number of spikes per trial and per electrode
    for i = 1:Nunits
        curr_unit = i;
        [spike_samples(:,i)] = squeeze(mean(get_spike_ts(curr_unit,stage,md)/get_dt));
    end
    
    spike_cs = spike_samples(correct_samples,:);
    
    
    if debug_on
        
        % Categorize samples the old fashioned way
        ctgsetli = get_ctgsetli(md);
        ctgsetli_cs = ctgsetli(correct_samples,:);
        ctgsetli_cs=ctgsetli_cs(:,1:8);
        
        % Estimate it instead based on conditions -> morphs -> ctgs
        ctgset_cs_est = conditions_to_trials(each_cs_condition);
        % Zero out the trials that lie along the 50% border regions
        index = sum(ctgset_cs_est,2);
        ctgset_cs_est(index == 1,:) = logical(ctgset_cs_est(index == 1,:)*0);

        % Compare category lists -- all good
        sum(sum(ctgsetli_cs - ctgset_cs_est))
        
        % Estimate firing rates for cells and categories
        Nctgs = size(ctgsetli_cs,2);
        for i = 1:Nctgs
            fr{i} = (spike_cs(ctgsetli_cs(:,i),:));
        end
        
        % Load old firing rate data for File 1
        sr = load('/Users/davestanley/Other/data/roy/buffer/analyze_data2/mode_5/stage_3/L011107.mat');
        
        % Plot firing rate comparison for single unit and ctg1
        figure; plot(sr.sfc{1}.numspikes(:,1))
        hold on; plot(fr{1}(:,1),'r')
        title('Unit1, Ctg1. ASSUMES WE ARE LOADING FILE #1, STAGE 3. If it doesnt match, likely wrong file'); legend('From Run analysis','Curr Estimate');
        xlabel('Trial in ctg1'); ylabel('Firing rate');
        
        % Compare all units from ctg1
        curr_ctg = 1; X = [mean(sr.sfc{curr_ctg}.numspikes); mean(fr{curr_ctg}); ]';
        figure; bar(X);
        title(['Mean firing rates for ctg ' num2str(curr_ctg)]);
        legend('From Run analysis','Curr Estimate');
        xlabel('Neuron #');
        diff(X,[],2)        % No difference
        
        % Compare morphs all cells to ctgs firing rates
        index = sch == 'A' & morphA > 55 & round(morphB) ~= 50; % Make sure to remove all border cases. Want to avoid rounding error
        frA = mean(spike_cs(index,:));
        index = sch == 'A' & morphA < 45 & round(morphB) ~= 50;
        frB = mean(spike_cs(index,:));
        figure; bar([frA ; mean(fr{1})]')
        frA - mean(fr{1})       % No difference
        
        
        
        
        % Plot firing rate vs morphs
        index = classify_scheme(sch,morphA,morphB,'A');
        curr_unit = 1;
        figure; plot(morphA(index),spike_cs(index,curr_unit),'.')
        x = morphA(index); y = spike_cs(index,curr_unit);
        [mean(y(x > 55)) mean(y(x < 45))]
        [mean(fr{1}(:,1)) mean(fr{2}(:,1))]
        
    end
    
    % Fixed the column ordering to match Jefferson's
    sens = zeros(Nunits,4);
    index = classify_scheme(sch,morphA,morphB,'A');
    
    if debug_on
        curr_neuron = 1;
        figure; scatter3(morphA(index)', morphB(index)', spike_cs(index,curr_neuron));
        xlabel('Morph A'); ylabel('Morph B'); zlabel('Firing Rate');
    end
    [sens(:,1)] = estimate_sensitivity(morphA(index),spike_cs(index,:));    % A relavent sensitivity
    [sens(:,4)] = estimate_sensitivity(morphB(index),spike_cs(index,:));    % A irrelevant sensitivity
    index = classify_scheme(sch,morphA,morphB,'B');
    [sens(:,2)] = estimate_sensitivity(morphB(index),spike_cs(index,:));    % B relavent sensitivity
    [sens(:,3)] = estimate_sensitivity(morphA(index),spike_cs(index,:));    % B irrelevant sensitivity
    
    
    unit_names = md.unit_names;
    unit_names=convert_unit_underscores(unit_names);
    
    % I GUESS I CAN LEAVE THE BAD TRIALS IN, SINCE THEY SHOULD NOT BE
    % AFFECTED BY THE CLIPPING

end

function ctgsetli = get_ctgsetli(md)


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
    
end



function [out] = estimate_sensitivity(morph,firing_rate)

    plot_on = 0;
    Nunits = size(firing_rate,2);
    
    x = morph(:);
    for i = 1:Nunits
        
        y = firing_rate(:,i);
        if plot_on; 
            Nbins = 60;
            uniques = unique(x);
            clear n bins
            for i = 1:length(uniques)
                [n(i,:), bins] = hist(y(x == uniques(i)),1:2:Nbins); 
            end
            figure; imagesc([uniques],[bins],n'); set(gca,'YDir','normal')
            xlabel('Morph percent'); ylabel('Firing rate'); title('Neuron firing rate histogram');
            colorbar;
        
        end
        if plot_on; figure; plot(x,y,'.'); xlabel('Morph %');ylabel('FiringRate');end
        
        p = polyfit(x,y,1);
        
        if plot_on; hold on; plot(x,polyval(p,x),'k','LineWidth',2); title(['Line slope m=' num2str(p(1))]); end
        out(i) = p(1);
        
    end
end



