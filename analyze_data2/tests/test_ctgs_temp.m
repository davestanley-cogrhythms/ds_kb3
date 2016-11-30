


% load('/Users/davestanley/Other/data/roy/CRC_Miller_data/O022805.mat')
load /Users/davestanley/Miller/roy/ctgX_merged_final/sansunits/O022805.mat
% load('/Users/davestanley/Other/data/roy/new_ctgX_files/O022805_new_conds_trials.mat')


Ntrials=length(start_times);

chosen_responses = [0 1 2 5 6];
chosen_responses = [0];

cri = logical(zeros(length(each_trial_response),1));

for i = 1:length(chosen_responses)
    cri = cri | (each_trial_response == chosen_responses(i));
end

cri = cri(:);
cri=cri';

good_ctg1 = build_logical(ctg1_trials,Ntrials) & cri
good_ctg2 = build_logical(ctg2_trials,Ntrials) & cri
good_ctg3 = build_logical(ctg3_trials,Ntrials) & cri
good_ctg4 = build_logical(ctg4_trials,Ntrials) & cri


% good_ctg1 = build_logical(ctg1_trials,Ntrials) & build_logical(correct_trials,Ntrials);
% good_ctg2 = build_logical(ctg2_trials,Ntrials) & build_logical(correct_trials,Ntrials);
% good_ctg3 = build_logical(ctg3_trials,Ntrials) & build_logical(correct_trials,Ntrials);
% good_ctg4 = build_logical(ctg4_trials,Ntrials) & build_logical(correct_trials,Ntrials);

% good_ctg1 = build_logical(ctg1_trials,Ntrials);
% good_ctg2 = build_logical(ctg2_trials,Ntrials);
% good_ctg3 = build_logical(ctg3_trials,Ntrials);
% good_ctg4 = build_logical(ctg4_trials,Ntrials);
% 
% n=1:Ntrials;
% y = ones(1,Ntrials);
% figure; plot(n,good_ctg1,'.');
% hold on; plot(n(ctg1_trials),y(ctg1_trials),'rx');
% 


ctgsetlis = [good_ctg1(:) good_ctg2(:) good_ctg3(:) good_ctg4(:)];


AA = sum(ctgsetlis)
AA(1)+AA(2)
AA(3)+AA(4)




