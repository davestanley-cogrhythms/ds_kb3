

% Testing script for test_filter60_lfp_sample_delay.m


%% Setup paths and load files
inpath = fullfile(getpath('path_lfp_sample_delay'));
inpath_temp = fullfile(inpath,'filtered_temp');
inpath_md = fullfile(inpath_temp,'filter_metadata');
lfp_files = get_filelist(inpath_temp);
metafiles = get_filelist(inpath_md);

%% Load and pack all metadata

md = cellfun(@(fn) load(fullfile(inpath_md,fn)), metafiles);

trials_per_file = arrayfun(@(s) size(s.power_reduction,1),md);
mintrials = min(trials_per_file);


% % Power reduction
clear s
s.d = {md.power_reduction};
s.mu = (cellfunu(@(d) mean(d), s.d)); s.mu = horzcat(s.mu{:});
s.std = (cellfunu(@(d) std(d), s.d)); s.std = horzcat(s.std{:});
s.all = (cellfunu(@(d) d(1:mintrials,:), s.d)); s.all = horzcat(s.all{:});
pr = s;


% % Interval lengths
clear s
s.d = {md.interval_lengths};
s.mu = (cellfunu(@(d) mean(d), s.d)); s.mu = horzcat(s.mu{:});
s.std = (cellfunu(@(d) std(d), s.d)); s.std = horzcat(s.std{:});
s.all = (cellfunu(@(d) d(1:mintrials,:), s.d)); s.all = horzcat(s.all{:});
il = s;

%% Load 60 Hz bad files (from identify_60Hz.m)
inpath = fullfile(getpath('path_badchannels'));
s = load(fullfile(inpath,'id_60Hz_s5.mat'));
bad_60hz = horzcat(s.bad_indices{:});
flnames = s.flnames;

% % Plotting
%% Plot 
n=1:length(pr.mu);

s = pr;
figure; subplot(121); hold on; plot(n,s.mu); %errorbar(n,s.mu,s.std); 
hold on; plot(n(bad_60hz),s.mu(bad_60hz),'r.'); legend('All files','Bad 60Hz'); ylabel('% Reduction post filter'); xlabel('Electrode number');
subplot(122); plot(n(~bad_60hz),s.mu(~bad_60hz)); legend('Good 60Hz files');

[fnum, unum] = unum2fnum2(695, flnames);


