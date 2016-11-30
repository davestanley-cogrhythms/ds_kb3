function [tsw tct] = get_switch_trials(md)

    Ntrials = length(md.start_times);

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
    %fprintf('Percent switch trials = %g \n', (Ntsw / (Ntsw + Ntct) *100) );
end