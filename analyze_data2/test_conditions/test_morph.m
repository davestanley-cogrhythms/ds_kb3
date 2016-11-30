% TEsts the performance of our condition2morph code

addpath('./funcs_supporting');

reload = 1;
remove_duplicate_conditions = 0; % I have assumed that each condition actually has 2 codes - 1 for match and 1 for non-match. Disregard one of these.

if reload
    load test_conds_f60.mat
end

etct = md.each_trial_condition;
% sample_on_trials = md.sample_on_trials;
Ntrials = length(etct);
ctg1t=build_logical(md.ctg1_trials,Ntrials);
ctg2t=build_logical(md.ctg2_trials,Ntrials);
ctg3t=build_logical(md.ctg3_trials,Ntrials);
ctg4t=build_logical(md.ctg4_trials,Ntrials);
ctgnr1t=build_logical(md.ctg1_nr_trials,Ntrials);
ctgnr2t=build_logical(md.ctg2_nr_trials,Ntrials);
ctgnr3t=build_logical(md.ctg3_nr_trials,Ntrials);
ctgnr4t=build_logical(md.ctg4_nr_trials,Ntrials);


ctgsett = [ctg1t' ctg2t' ctg3t' ctg4t' ctgnr1t' ctgnr2t' ctgnr3t' ctgnr4t'];

ctgsett_est = conditions_to_trials(etct);

% Zero out the trials that lie along the 50% border regions
index = sum(ctgsett_est,2);
ctgsett_est(index == 1,:) = logical(ctgsett_est(index == 1,:)*0);

% Plot overlay of ctgs vs ctg_ests. Our morph code is working!
figure; imagesc(ctgsett);
hold on;
index = find(ctgsett_est == 1);
[I,J] = ind2sub(size(ctgsett_est),index);
hold on; plot(J,I,'wx');



