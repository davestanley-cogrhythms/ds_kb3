
%% Load a file

clear
close all

addpath('./funcs_supporting')

path_md = getpath('path_metadata');
file_list = get_filelist(path_md);

fnum = 49;
fnum = 5;
fnum = 10;

md = load(fullfile(path_md,file_list{fnum}));
vars_pull(md);


%% Define constants
Npoints = 16;
roys_way = 0;       % Roy's way is estimating based on conditions being even/odd. My
                    % way uses all encodes. My way is better for Cromer
                    % data because there are a few trials that for some
                    % reason aren't considered error trials, even though
                    % they have some weird all_encodes patterns (and
                    % severely long response times)

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

% Alias for time
t = all_encode_times;


%% Some common variables

sample_on_encodes=find(all_encodes == code_sample_on);   % Sample on encode indices

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




%% Break data up into short and long trials (Roy's way)


ae_inds = sample_on_encodes(inds);
ae_indl = sample_on_encodes(indl);
ae_ind_rest = sample_on_encodes(ind_rest);
ae_ind_err = sample_on_encodes(ind_err);

figure;


[indg] = get_centered_encodes(ae_inds,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,ae_inds,all_encodes,t);
hold on; h1=plot(tcent,aesq,'b');

[indg] = get_centered_encodes(ae_indl,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,ae_indl,all_encodes,t);
hold on; h2=plot(tcent,aesq,'r');

[indg] = get_centered_encodes(ae_ind_rest,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,ae_ind_rest,all_encodes,t);
hold on; h3=plot(tcent,aesq,'k');

[indg] = get_centered_encodes(ae_ind_err,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,ae_ind_err,all_encodes,t);
hold on; h4=plot(tcent,aesq-100,'g');

legend([h1(1) h2(1) h4(1)], 'Match','Non-match','Error')
ylabel('All encodes'); xlabel('Time (ms) centered around sample on');

%% Plot encodes time vs diff time
    % For Roy, most variability happens 3 data point after sample on for
    % match; 5 for non-match
    % For Cromer, it's actually the same. The reason Cromer's code shift is
    % chosen to be 4 is because there are a few "weird" trials in cromer
    % that I want to avoid. See my email to Jefferson (search BU email for
    % keyword all_encodes or for subject match/non-match)
    
    
figure;
[indg] = get_centered_encodes(ae_inds,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,ae_inds,all_encodes,t);
hold on; h1=plot((1:Npoints-1)'-7,diff(tcent,[],1),'b');

[indg] = get_centered_encodes(ae_indl,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,ae_indl,all_encodes,t);
hold on; h2=plot((1:Npoints-1)'-7,diff(tcent,[],1),'r');

[indg] = get_centered_encodes(ae_ind_rest,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,ae_ind_rest,all_encodes,t);
if ~isempty(tcent); hold on; h3=plot((1:Npoints-1)'-7,diff(tcent,[],1),'k'); end

% [indg] = get_centered_encodes(ind_err,Npoints);
% [aesq, tcent] = lookup_plotting_data(indg,ind_err,all_encodes,t);
% hold on; h4=plot((1:Npoints-1)'-7,diff(tcent,[],1),'g');

legend([h1(1) h2(1)], 'Match','Non-match')
ylabel('diff(all encode times)'); xlabel('All encode times index (allencodes sample on = 0)');





%% Plot reaction times


%match_dt = all_encode_times(ind2+code_match_response_shift-1)-all_encode_times(ind2); % test time since sample on - should be about 1.6 seconds
match_rt = all_encode_times(sample_on_encodes+code_match_response_shift)-all_encode_times(sample_on_encodes+code_match_response_shift-1);
match_rt(~inds) = NaN;

%nonmatch_dt = all_encode_times(ind2+code_nonmatch_response_shift-1)-all_encode_times(ind2); % Test time since sample on - should be about 3.2s
nonmatch_rt = all_encode_times(sample_on_encodes+code_nonmatch_response_shift)-all_encode_times(sample_on_encodes+code_nonmatch_response_shift-1);
nonmatch_rt(~indl) = NaN;

figure; hist(match_rt,20); title('Match RT');
figure; hist(nonmatch_rt,20); title('Nonmatch RT');


