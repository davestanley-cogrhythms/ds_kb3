
function samplesli = get_good_samples(trials,md)

    plot_debug = 0;
    exclude_incorrect_trials = 1;
    corrtrcd = 0;       % "Code" for correct trial. Recall, 0 = Monkey was correct.

    %vars_pull(md);      % Pull variables out of structure and place in local work space - this command is slow, so don't use it!
    each_trial_response = md.each_trial_response;
    sample_on_trials = md.sample_on_trials;
    
    Ntrials = length(md.start_times);  % Total number of trials
    
    %good_ind = sum(isnan(C_arr(:,:,unit_to_plot)),1) == 0;
    chosen_ctg_trials = build_logical(trials,Ntrials);
    
    if exclude_incorrect_trials
        chosen_trials = chosen_ctg_trials(:) & each_trial_response(:) == corrtrcd;  % And chosen category with "correct" trials
    else
        chosen_trials = chosen_ctg_trials(:) & sample_on_trials(:);  % And chosen category with "correct" trials
    end
    
    
    samplesli = (trials_to_samples( chosen_trials, sample_on_trials));  % Convert chosen trials from being in terms of trials to being in terms of sample ons
    
    if plot_debug
        figure;
        n=1:length(sample_on_trials);
        plot(n,chosen_ctg_trials(:) & each_trial_response(:) == corrtrcd,'b.');
        hold on;
        n_sample_ons = n(sample_on_trials);
        plot(n_sample_ons,samplesli,'gx');
        
        x = n_sample_ons(samplesli);
        y = n(chosen_trials);
        
        sum(x ~= y)
        
        legend('Chosen Trials in Trial Space','Chosen Trials in Sample On Space');
        
    end
end

