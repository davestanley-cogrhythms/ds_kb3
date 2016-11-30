
%% Reaction times all files
why
why
why

clear
% close all

addpath(genpath('../../funcs_supporting_local/'));

path_md = fullfile(getpath('path_metadata'),'sansunits');
%path_md = fullfile(getpath('path_metadata'));
%path_md = fullfile(getpath('path_roydata_orig'),'metadata');
%path_md = getpath('path_roy_newctgX');

file_list = get_filelist(path_md);

md = cellfun(@load,fullfile(path_md,file_list));


%% Define constants

roys_way = 0;

if ~get_iscromer
    code_shift = 3;
    code_match_response_shift = 3;
    code_nonmatch_response_shift = 5; % I think
    code_sample_on = 27;
    code_test = 29;
    code_short = 4;
    code_long = 30;
else
%     % Old method ... this one catches some bugs
%     code_shift = 3;
%     code_sample_on = 23;
%     code_test = 25;
%     code_short = 175;
%     code_long = 26;
    
    code_shift = 4;
    code_match_response_shift = 3;
    code_nonmatch_response_shift = 5;
    code_sample_on = 23;
    code_test = 25;
    code_short = 176;
    code_long = 27;
end



%% Run analysis

clear match_rt nonmatch_rt schA schB schA2
for i = 1:length(md)
    vars_pull(md(i));
    
    
    % Find all instances of 27 - these are sample on indices
    sample_on_encodes=find(all_encodes == code_sample_on);


    %Error trials
    etrs = each_trial_response(sample_on_trials);
    ind_err = etrs ~= 0;     % Trials for which the correct answer was not given (for whatever reason)

    if roys_way
        inds = logical(mod(each_trial_condition,2));        % Odd for match
        indl = ~logical(mod(each_trial_condition,2));       % Even for non-match
        inds = inds(sample_on_trials) & ~ind_err;
        indl = indl(sample_on_trials) & ~ind_err;
    else % My way, using all_encodes
        encodesp4 = all_encodes(sample_on_encodes+code_shift);
        inds = encodesp4 == code_short & ~ind_err;
        indl = encodesp4 == code_long & ~ind_err;
    end

    ind_rest = ~(inds | indl | ind_err);

    %match_dt = all_encode_times(ind2+code_match_response_shift-1)-all_encode_times(ind2); % test time since sample on - should be about 1.6 seconds
    match_rt_calc{i} = all_encode_times(sample_on_encodes+code_match_response_shift)-all_encode_times(sample_on_encodes+code_match_response_shift-1);
    match_rt_calc{i}(~inds | ind_err) = NaN;

    %nonmatch_dt = all_encode_times(ind2+code_nonmatch_response_shift-1)-all_encode_times(ind2); % Test time since sample on - should be about 3.2s
    nonmatch_rt_calc{i} = all_encode_times(sample_on_encodes+code_nonmatch_response_shift)-all_encode_times(sample_on_encodes+code_nonmatch_response_shift-1);
    nonmatch_rt_calc{i}(~indl | ind_err) = NaN;
    
    % Get current schemes
    Ntrials = length(start_times);
%         % Slow way
%     clear sch morphA morphB
%     for j = 1:Ntrials
%         [sch(j),morphA(j),morphB(j)] = condition2morph(each_trial_condition(j)); % Slow, but not limiting.
%     end
%     schA2{i} = get_good_samples(find(sch == 'A' & morphA ~= 50 & morphB ~= 50),md(i));
    
        % Faster way
    schA_tr = [ctg1_trials; ctg2_trials];
    schA{i} = get_good_samples(schA_tr,md(i))';
    schB_tr = [ctg3_trials; ctg4_trials];
    schB{i} = get_good_samples(schB_tr,md(i))';
    
end


%% Plot all
figure
for i = 1:length(md)
    clf;
    subplot(211); hist(match_rt_calc{i}); title('Match');        % Match
    subplot(212); hist(nonmatch_rt_calc{i}); title('Non-match');     % Non-match
    pause
end


%% Calculate RT stats for all trials, schA, and schB

stats_fun = @mean;
% stats_fun = @median;

if ~get_iscromer
    range_all=1:79;
    range_monkey1=1:40;
    range_monkey2=41:79;
else
    temp={md.recording_name};
    range_all=1:length(md);
    range_monkey1=find(cellfun(@(x) ~isempty(x),strfind(temp,'lu')));
    range_monkey2=find(cellfun(@(x) ~isempty(x),strfind(temp,'ti')));
end

medRT = cellfun(@(x) stats_fun(x(~isnan(x))),match_rt_calc);
medRT_A = cellfun(@(x,y) stats_fun(x(~isnan(x) & y)),match_rt_calc,schA);    % Mean reaction time for trials that are (1) match trials (i.e. non NaNs) and (2) Scheme A
medRT_B = cellfun(@(x,y) stats_fun(x(~isnan(x) & y)),match_rt_calc,schB);


%% Plot Monkey L vs Monkey O, all trials

figure;
data1 = medRT(range_monkey1);
data2 = medRT(range_monkey2);

% Plot data1 and data2
hold on; plot(data1,'b')
hold on; plot(data2,'r')

% Stats data1 and data2
med1 = stats_fun(data1);
med2 = stats_fun(data2);
[p, h] = ranksum(data1,data2);
% [p, h] = signrank(data1,data2);

legend(['Monkey L mean =' num2str(med1)],['Monkey O mean =' num2str(med2)])
xlabel('File'); ylabel('Med Reaction Time')
title(['Monkey L vs O RT; ranksum p=' num2str(p)])




%% Plot Monkey all AvB

figure;
data1 = medRT_A(range_all);
data2 = medRT_B(range_all);

% Plot data1 and data2
hold on; plot(data1,'b')
hold on; plot(data2,'r')

% Stats data1 and data2
med1 = stats_fun(data1);
med2 = stats_fun(data2);
[p, h] = ranksum(data1,data2);
[p, h] = signrank(data1,data2);

legend(['Sch A mean =' num2str(med1)],['Sch B mean =' num2str(med2)])
xlabel('File'); ylabel('Med Reaction Time')
title(['Monkey All RT; signrank p=' num2str(p)])

%% Plot Monkey L AvB

figure;
data1 = medRT_A(range_monkey1);
data2 = medRT_B(range_monkey1);

% Plot data1 and data2
hold on; plot(data1,'b')
hold on; plot(data2,'r')

% Stats data1 and data2
med1 = stats_fun(data1);
med2 = stats_fun(data2);
[p, h] = ranksum(data1,data2);
[p, h] = signrank(data1,data2);

legend(['Sch A mean =' num2str(med1)],['Sch B mean =' num2str(med2)])
xlabel('File'); ylabel('Med Reaction Time')
title(['Monkey L RT; signrank p=' num2str(p)])


%% Plot Monkey O AvB

figure;
data1 = medRT_A(range_monkey2);
data2 = medRT_B(range_monkey2);

% Plot data1 and data2
hold on; plot(data1,'b')
hold on; plot(data2,'r')

% Stats data1 and data2
med1 = stats_fun(data1);
med2 = stats_fun(data2);
[p, h] = ranksum(data1,data2);
[p, h] = signrank(data1,data2);

legend(['Sch A mean =' num2str(med1)],['Sch B mean =' num2str(med2)])
xlabel('File'); ylabel('Med Reaction Time')
title(['Monkey O RT; signrank p=' num2str(p)])


