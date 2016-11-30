

% This code is used to verify that the convert_to_seconds part of
% extract_metadata_merge is working correctly.


clear
close all


path_md = fullfile(getpath('path_metadata'),'sansunits');
%path_md = fullfile(getpath('path_metadata'));
%path_md = fullfile(getpath('path_roydata_orig'),'metadata');
%path_md = getpath('path_roy_newctgX');

file_list = get_filelist(path_md);

md = cellfun(@load,fullfile(path_md,file_list));


%% Plot some spacing between sample ons & and all_encode_times
[aet,sot,etr,sotr] = struct_pull(md,{'all_encode_times','sample_on_times','each_trial_response','sample_on_trials'});

etr = cellfunu(@(etr,sotr) etr(sotr),etr,sotr);
etr = cellfunu(@(etr) etr(1:end-1),etr);

dsot = cellfunu(@(x) (diff(x)),sot);
dsot_good = cellfunu(@(dsot,etr) dsot(etr == 0),dsot,etr);  % Diffs in sample on times only for trials that were correct (error code 0)

mu = cellfun(@(dsot) median(dsot),dsot);
figure; plot(mu)
xlabel('File Number'); ylabel('Median trial duration');

mu = cellfun(@(dsot) mean(dsot),dsot_good);
figure; plot(mu)
xlabel('File Number'); ylabel('Mean good trial duration');


