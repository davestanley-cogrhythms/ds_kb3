
%save(fullfile(getpath('path_matfiles_store'),'md_temp.md'),'md')

clear
% close all

addpath('./funcs_supporting')

path_md = getpath('path_metadata');
file_list = get_filelist(path_md);

fnum = 49;
fnum = 5;

md = load(fullfile(path_md,file_list{fnum}));
vars_pull(md);

separatrix = 9;

%% Guess all_encodes code for sample on
sample_on_times=sample_on_times*1000;
Ntr = length(sample_on_times);
ind = zeros(1,Ntr);
for i = 1:Ntr
    ind(i) = find(all_encode_times == sample_on_times(i)/1e3);
end

% The all encodes correspondign to sample on times
sample_codes = all_encodes(ind);
figure; hist(sample_codes)
title('Sample on code is 27 for Roy, 23 for Cromer');


%% Define constants
if ~get_iscromer
    code_shift = 3;
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
    code_sample_on = 23;
    code_test = 25;
    code_short = 176;
    code_long = 27;
end

% Find all instances of 27 - these are sample on indices
sample_on_encodes=find(all_encodes == code_sample_on);
N = length(all_encodes);
n=1:N;
t = all_encode_times;



%% Plot all instances of 27
figure;
plot(n,all_encodes,'k');
hold on; plot(n(sample_on_encodes),all_encodes(sample_on_encodes),'kx');
xlabel('All encodes Index'); ylabel('All encodes'); legend('All encodes','Sample On Code (27)')


ind3=find(all_encodes == code_short);
hold on; plot(n(ind3),all_encodes(ind3),'b.');
ind3=find(all_encodes == code_long);
hold on; plot(n(ind3),all_encodes(ind3),'r.');
legend('encodes','sample','code 4','code 30')

%% Plot all encodes versus mod all encodes - collapsed into one trial; assume each all encodes is length 18
% This not working right
%figure; plot(n,all_encodes);
% clf; hold on; plot(mod(n(1:30),24),all_encodes(1:30),'r');



%% Now plot all_encodes and overlay with sample_on_times and test_on_times
figure;
plot(all_encode_times,all_encodes);
ab = ones(1,Ntr);
hold on; plot(sample_on_times/1e3,code_sample_on*ab,'gx');
hold on; plot(sample_on_times/1e3+1.6,code_test*ab,'rx');
hold on; plot(all_encode_times(sample_on_encodes+2),all_encodes(sample_on_encodes+2),'mx');
xlabel('All encode times'); ylabel('All encodes'); legend('All encodes','Val of 27','Estimated test on','2 events after sample on')
% Seems to be about right - test image is the code 2 after sample on

%% Scatterplot time difference versus encodes value; also hist of all encodes values
figure;
subplot(211);plot(all_encode_times(sample_on_encodes+2) - all_encode_times(sample_on_encodes), all_encodes(sample_on_encodes+2),'.'); xlabel('Time difference'); ylabel('All Encodes');
subplot(212);hist(all_encodes(sample_on_encodes+2)); xlabel('Encodes value');
% This index does not always appear the same time after sample on, nor
% does it always equal 29; however, this could be due to failed trials.

%% Distribution of trial lengths
%figure; hist(diff(ind2),1000)
figure; hist(diff(all_encode_times(sample_on_encodes)),1000); xlabel('Sample on time differences.'); %xlim([7e3 11e3]);

%% Plot chunks of data centered around sample on

Npoints = 16;

% Error trials
etrs = each_trial_response(sample_on_trials);
ind_err = etrs ~= 0;     % Trials for which the correct answer was not given (for whatever reason)
ind_good = sample_on_encodes(~ind_err);

[indg] = get_centered_encodes(ind_good,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,ind_good,all_encodes,t);
figure; plot(tcent,aesq);


%% Break data up into short and long trials (Roy's way)

Npoints = 16;

%Error trials
etrs = each_trial_response(sample_on_trials);
ind_err = etrs ~= 0;     % Trials for which the correct answer was not given (for whatever reason)

inds = logical(mod(each_trial_condition,2));        % Odd for match
indl = ~logical(mod(each_trial_condition,2));       % Even for non-match

inds = inds(sample_on_trials) & ~ind_err;
indl = indl(sample_on_trials) & ~ind_err;

ind_rest = ~(inds | indl | ind_err);

inds = sample_on_encodes(inds);
indl = sample_on_encodes(indl);
ind_rest = sample_on_encodes(ind_rest);
ind_err = sample_on_encodes(ind_err);

figure;


[indg] = get_centered_encodes(inds,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,inds,all_encodes,t);
hold on; h1=plot(tcent,aesq,'b');

[indg] = get_centered_encodes(indl,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,indl,all_encodes,t);
hold on; h2=plot(tcent,aesq,'r');

[indg] = get_centered_encodes(ind_rest,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,ind_rest,all_encodes,t);
hold on; h3=plot(tcent,aesq,'k');

[indg] = get_centered_encodes(ind_err,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,ind_err,all_encodes,t);
hold on; h4=plot(tcent,aesq-100,'g');

legend([h1(1) h2(1) h4(1)], 'Match','Non-match','Error')
ylabel('All encodes'); xlabel('Time (ms) centered around sample on');



%% Break data up into short and long trials (new way)
% Find all short and long sample on trials; assume these are
% match/non-match

%Error trials
etrs = each_trial_response(sample_on_trials);
ind_err = etrs ~= 0;     % Trials for which the correct answer was not given (for whatever reason)

encodesp4 = all_encodes(sample_on_encodes+code_shift);
inds = encodesp4 == code_short & ~ind_err;
indl = encodesp4 == code_long & ~ind_err;
ind_rest = ~(inds | indl | ind_err);


etcs = each_trial_condition(sample_on_trials);
figure; plot(mod(etcs(inds),2));
hold on; plot(mod(etcs(indl),2),'r');
xlabel('trial #'); ylabel('even vs odd each trial condition'); legend('short trials','long trials');

cright= mod(etcs(inds),2) == 1;
cwrong = mod(etcs(inds),2) == 0;

n = 1:length(etcs);
figure; plot(n,etcs);
etcss = etcs(inds);
ns = n(inds);
hold on; plot(ns(cright),etcss(cright),'g.');
hold on; plot(ns(cwrong),etcss(cwrong),'r.');
xlabel('trial #'); ylabel('each trial condition'); legend('All shorts','Shorts are odd','Shorts even (classifier fail)');
successful_trial_conditions = etcss(cright)'
failed_trial_conditions = etcss(cwrong)'

% sstca = [sstca successful_trial_conditions];
% sftca = [sftca failed_trial_conditions];
% 
% sstca = sort(unique(sstca));
% sftca = sort(unique(sftca));
% 
% arrayfun(@(x) sum(x == sstca),sftca)
% 
% etcsl = etcs(mod(etcs,2) == 0);
% etcsl = sort(unique(etcsl));
% arrayfun(@(x) sum(x == etcsl),sftca)



inds = sample_on_encodes(inds);
indl = sample_on_encodes(indl);
ind_rest = sample_on_encodes(ind_rest);
ind_err = sample_on_encodes(ind_err);

figure; 

[indg] = get_centered_encodes(inds,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,inds,all_encodes,t);
hold on; h1=plot(tcent,aesq,'b');

[indg] = get_centered_encodes(indl,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,indl,all_encodes,t);
hold on; h2=plot(tcent,aesq,'r');

[indg] = get_centered_encodes(ind_rest,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,ind_rest,all_encodes,t);
hold on; h3=plot(tcent,aesq,'k');

[indg] = get_centered_encodes(ind_err,Npoints);
[aesq, tcent] = lookup_plotting_data(indg,ind_err,all_encodes,t);
hold on; h4=plot(tcent,aesq-100,'g');

legend([h1(1) h2(1) h4(1)], 'Match','Non-match','Error')
ylabel('All encodes'); xlabel('Time (ms) centered around sample on');


%% Break data up into short and long trials (old way)
% % Find all short and long sample on trials; assume these are
% % match/non-match
% 
% %Error trials
% etrs = each_trial_response(sample_on_trials);
% ind_err = etrs ~= 0;     % Trials for which the correct answer was not given (for whatever reason)
% 
% dsot = diff(all_encode_times(ind2));  %Diff in sample on times (also equals diff(sample_on_times))
% inds = (dsot < separatrix); inds = [inds; false]; inds = inds & ~ind_err;     % Remove error trials
% indl = (dsot >= separatrix & dsot < separatrix+2); indl = [indl;false]; indl = indl & ~ind_err;
% ind_rest = ~(inds | indl | ind_err);
% inds = ind2(inds);
% indl = ind2(indl);
% ind_rest = ind2(ind_rest);
% ind_err = ind2(ind_err);
% 
% [indg] = get_centered_encodes(inds,Npoints);
% [aesq, tcent] = lookup_plotting_data(indg,inds,all_encodes,t);
% figure; h1=plot(tcent,aesq,'b');
% 
% 
% [indg] = get_centered_encodes(indl,Npoints);
% [aesq, tcent] = lookup_plotting_data(indg,indl,all_encodes,t);
% hold on; h2=plot(tcent,aesq,'r');
% 
% 
% [indg] = get_centered_encodes(ind_rest,Npoints);
% [aesq, tcent] = lookup_plotting_data(indg,ind_rest,all_encodes,t);
% hold on; h3=plot(tcent,aesq,'k');
% 
% [indg] = get_centered_encodes(ind_err,Npoints);
% [aesq, tcent] = lookup_plotting_data(indg,ind_err,all_encodes,t);
% hold on; h4=plot(tcent,aesq-100,'g');
% 
% legend([h1(1) h2(1) h3(1) h4(1)], 'Match','Non-match','Rest','Error')
% ylabel('All encodes'); xlabel('Time (ms) centered around sample on');



