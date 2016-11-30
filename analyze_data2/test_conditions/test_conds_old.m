
reload = 1;

if reload
    load test_conds_f1.mat
    ctgsetli = ctgsetli(:,1:8);
end

each_trial_condition = md.each_trial_condition;
sample_on_trials = md.sample_on_trials;


etc_samp = each_trial_condition(sample_on_trials);

%Scheme A
sch_cat{1} = ctgsetli(:,1) & ctgsetli(:,7);% A_AB11
sch_cat{2}= ctgsetli(:,1) & ctgsetli(:,8);% A_AB12
sch_cat{3} = ctgsetli(:,2) & ctgsetli(:,7);% A_AB21
sch_cat{4} = ctgsetli(:,2) & ctgsetli(:,8);% A_AB22

%Scheme B
sch_cat{5} = ctgsetli(:,3) & ctgsetli(:,5);% B_AB11
sch_cat{6} = ctgsetli(:,4) & ctgsetli(:,5);% B_AB12
sch_cat{7} = ctgsetli(:,3) & ctgsetli(:,6);% B_AB21
sch_cat{8} = ctgsetli(:,4) & ctgsetli(:,6);% B_AB22

figure; imagesc(etc_samp);
N=length(sch_cat);
for i = 1:N
    codes_in_state{i} = etc_samp(sch_cat{i});
    unique_codes{i} = unique(codes_in_state{i});
    subplotsq(N,i); hist(codes_in_state{i});
    unique_codes{i}'
end

all_unique = sort(cat(1,unique_codes{:}));
all_unique2 = unique(etc_samp);
all_unique3 = unique(each_trial_condition);

all_cat = any(cat(2,sch_cat{:}),2);
no_cat = ~all_cat;

codes_nostate = etc_samp(no_cat);
unique_nostate_codes = unique(codes_nostate);

index = find(etc_samp == 13); ctgsetli(index,:)
index  = find(each_trial_condition == 14)



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





