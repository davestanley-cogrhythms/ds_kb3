
reload = 1;
remove_duplicate_conditions = 0; % I have assumed that each condition actually has 2 codes - 1 for match and 1 for non-match. Disregard one of these.

if reload
    load test_conds_f1.mat
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

%Scheme A
sch_cat{1} = ctgsett(:,1) & ctgsett(:,7);% A_AB11
sch_cat{2}= ctgsett(:,1) & ctgsett(:,8);% A_AB12
sch_cat{3} = ctgsett(:,2) & ctgsett(:,7);% A_AB21
sch_cat{4} = ctgsett(:,2) & ctgsett(:,8);% A_AB22

%Scheme B
sch_cat{5} = ctgsett(:,3) & ctgsett(:,5);% B_AB11
sch_cat{6} = ctgsett(:,4) & ctgsett(:,5);% B_AB12
sch_cat{7} = ctgsett(:,3) & ctgsett(:,6);% B_AB21
sch_cat{8} = ctgsett(:,4) & ctgsett(:,6);% B_AB22

figure; imagesc(ctgsett);
N=length(sch_cat);
for i = 1:N
    codes_in_state{i} = etct(sch_cat{i});
    unique_codes{i} = unique(codes_in_state{i});
    subplotsq(N,i); hist(codes_in_state{i});
    if remove_duplicate_conditions
        unique_codes{i} = unique_codes{i}(mod(unique_codes{i},2) == 0) / 2;     % Reduce pairs of entries to single entries (removes the match/notmatch duplication)
    end
    sort(unique_codes{i})'
    
end


% % Base unique codes off of md.ctg instad of trials
clear unique_codes
for i = 1:4
    eval(['unique_codes{i} = sort(md.ctg' num2str(i) '_conds);']);
    %unique_codes{i}=unique_codes{i}';
    if remove_duplicate_conditions
        unique_codes{i} = unique_codes{i}(mod(unique_codes{i},2) == 0) / 2;     % Reduce pairs of entries to single entries (removes the match/notmatch duplication)
    end
    sort(unique_codes{i})'
    
end

% Plot unique codes in trials with assigned ctgs vs unique codes in all trials
all_unique = sort(cat(1,unique_codes{:}))'
all_unique2 = unique(etct)'


% Look at all codes that do NOT have an assigned ctg
all_cat = any(cat(2,sch_cat{:}),2);
no_cat = ~all_cat;

codes_nostate = etct(no_cat);
unique_nostate_codes = unique(codes_nostate)'
if remove_duplicate_conditions
    unique_nostate_codes = unique_nostate_codes(mod(unique_nostate_codes,2) == 0) / 2
end

index = find(etct == 13); ctgsett(index,:)
index  = find(etct == 14)


% Verify unique codes - these should be equal
unique(etct(ctg1t))' == md.ctg1_conds'



% 
% %SCheme A
% sch_border{1} = ctgsetli(:,1) & ~ctgsetli(:,7) & ~ctgsetli(:,8);
% sch_border{2} = ctgsetli(:,2) & ~ctgsetli(:,7) & ~ctgsetli(:,8);
% 
% % Scheme B
% sch_border{3} = ctgsetli(:,3) & ~ctgsetli(:,5) & ~ctgsetli(:,6);
% sch_border{4} = ctgsetli(:,4) & ~ctgsetli(:,5) & ~ctgsetli(:,6);
% 
% 
% figure;
% N=length(sch_border);
% for i = 1:N
%     codes_in_state{i} = etc_samp(sch_border{i});
%     unique_codes{i} = unique(codes_in_state{i});
%     subplotsq(N,i); hist(codes_in_state{i});
%     unique_codes{i}'
% end





