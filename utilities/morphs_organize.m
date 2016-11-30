
function morphs_organize
    % Organizes files so that the .con files match to the 1:79 data files

    % Setup paths
    path_roy = getpath('path_roy');
    outpath = getpath('path_conditions');

    % Load congruent/non-congruent trials
    load(fullfile(getpath('path_conditions'),'conditions.mat'));
    
    % Get expected file list
    path_metadata = fullfile(getpath('path_metadata'),'sansunits');
    file_list = get_filelist(path_metadata);
    
    ind = zeros(1,length(file_list));
    
    for i = 1:length(file_list)
        ind(i) = filename2cond_index(fnames,file_list{i});
    end
    fnames = fnames(ind);
    textdata=textdata(ind);
    conds_to_test1 = conds_to_test1(ind);
    conds_to_test2 = conds_to_test2(ind);
    
    
    
    % Set up data for testing
    md_all = cellfun(@load,fullfile(path_metadata,file_list));
    
    for i = 1:length(file_list);
        %%
        fnum = i;
        md = md_all(fnum);

        condsmap = conds_to_test1{fnum};
        
        iscat = condsmap(md.each_trial_condition,2);
        isthin = condsmap(md.each_trial_condition,3);
        
        if ~get_iscromer
            if fnum > 40     % For monkey O, numbering is reversed
                isthin = ~isthin;
            end
        end
        
        % Make sure that fat/thin classificaiton matches up with existing metadata; return a warning if it fails
        Ntrials = length(md.start_times);  % Total number of trials
        ctg1_trials = build_logical(md.ctg1_trials,Ntrials);
        ctg2_trials = build_logical(md.ctg2_trials,Ntrials);
        ctg3_trials = build_logical(md.ctg3_trials,Ntrials);
        ctg4_trials = build_logical(md.ctg4_trials,Ntrials);
        ctg1_test = (ctg1_trials & md.tr_match') | (ctg2_trials & md.tr_nonmatch');
        ctg2_test = (ctg2_trials & md.tr_match') | (ctg1_trials & md.tr_nonmatch');
        ctg3_test = (ctg3_trials & md.tr_match') | (ctg4_trials & md.tr_nonmatch');
        ctg4_test = (ctg4_trials & md.tr_match') | (ctg3_trials & md.tr_nonmatch');
        if sum(ctg1_test(:)) ~= sum(ctg1_test(:) & iscat(:)); warning('Incorrect classification of train images');end
        if sum(ctg2_test(:)) ~= sum(ctg2_test(:) & ~iscat(:)); warning('Incorrect classification of train images');end
        if sum(ctg3_test(:)) ~= sum(ctg3_test(:) & isthin(:)); warning('Incorrect classification of train images');end
        if sum(ctg4_test(:)) ~= sum(ctg4_test(:) & ~isthin(:)); warning('Incorrect classification of train images');end
        
        ctg_all = get_good_samples(1:length(md.sample_on_trials),md);   % I am updating ctg_all to only include correct trials. Prior to analyze_data7, this also included incorrect trials.
        ctg_all = ctg_all(:);
    

        % Find congruent vs incongruent trials
        ctg1_nr_trials = build_logical(md.ctg1_nr_trials,Ntrials);
        ctg2_nr_trials = build_logical(md.ctg2_nr_trials,Ntrials);
        ctg3_nr_trials = build_logical(md.ctg3_nr_trials,Ntrials);
        ctg4_nr_trials = build_logical(md.ctg4_nr_trials,Ntrials);
        
        % For match case
            % Congruent: testing images match
        match_congruent{i} = (md.tr_match(:) & ctg1_nr_trials(:) & iscat(:)) | ... 
                            (md.tr_match(:) & ctg2_nr_trials(:) & ~iscat(:)) | ... ;
                            (md.tr_match(:) & ctg3_nr_trials(:) & isthin(:)) | ... ;
                            (md.tr_match(:) & ctg4_nr_trials(:) & ~isthin(:)) ;
            % Incongruent: testing images non-match
        match_incongruent{i} = (md.tr_match(:) & ctg1_nr_trials(:) & ~iscat(:)) | ... 
                            (md.tr_match(:) & ctg2_nr_trials(:) & iscat(:)) | ... ;
                            (md.tr_match(:) & ctg3_nr_trials(:) & ~isthin(:)) | ... ;
                            (md.tr_match(:) & ctg4_nr_trials(:) & isthin(:)) ;
        
        % For nonmatch case
            % Congruent: testing images non-match
        nonmatch_congruent{i} = (md.tr_nonmatch(:) & ctg1_nr_trials(:) & ~iscat(:)) | ... 
                            (md.tr_nonmatch(:) & ctg2_nr_trials(:) & iscat(:)) | ... ;
                            (md.tr_nonmatch(:) & ctg3_nr_trials(:) & ~isthin(:)) | ... ;
                            (md.tr_nonmatch(:) & ctg4_nr_trials(:) & isthin(:)) ;
            % Incongruent: testing images non-match
        nonmatch_incongruent{i} = (md.tr_nonmatch(:) & ctg1_nr_trials(:) & iscat(:)) | ... 
                            (md.tr_nonmatch(:) & ctg2_nr_trials(:) & ~iscat(:)) | ... ;
                            (md.tr_nonmatch(:) & ctg3_nr_trials(:) & isthin(:)) | ... ;
                            (md.tr_nonmatch(:) & ctg4_nr_trials(:) & ~isthin(:)) ;
        
        
    end
    
    save(fullfile(outpath,'conditions_organized.mat'),'fnames','conds_to_test1','conds_to_test2','textdata','info_str','match_congruent','match_incongruent','nonmatch_congruent','nonmatch_incongruent');
    
end
