

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
[aet,sot] = struct_pull(md,{'all_encode_times','sample_on_times'});




mus = cellfun(@(x) median(diff(x)),sot);
figure; plot(mus)


mus = cellfun(@(x) median(diff(x)),aet);
figure; plot(mus/1e0)

%% Plot some spacing between units (note - need to use the non-sans units version for this)
if isempty(strfind(path_md,'sansunits'))    % This code only works when we load the non-sansunits md files
    [unit_array] = struct_pull(md,{'unit_array'});

    unit_array = horzcat(unit_array{:});
    hist(diff(unit_array{1}),100);
    mus = cellfun(@(x) max((x)),unit_array);
    figure; plot(mus)

end


%% Plot all the switches indicating whether stats were converted

[soc] = struct_pull(md,'metadata');
[soc] = struct_pull(soc,'units_converted');
[soc] = struct_pull(soc,'sample_on');
soc = cell2mat(soc);
figure; plot(soc)

[soc] = struct_pull(md,'metadata');
[soc] = struct_pull(soc,'units_converted');
[soc] = struct_pull(soc,'unit_array');
soc = cell2mat(soc);
figure; plot(soc)

[soc] = struct_pull(md,'metadata');
[soc] = struct_pull(soc,'units_converted');
[soc] = struct_pull(soc,'start_stop_times');
soc = cell2mat(soc);
figure; plot(soc);


[soc] = struct_pull(md,'metadata');
[soc] = struct_pull(soc,'units_converted');
[soc] = struct_pull(soc,'all_encode_times');
soc = cell2mat(soc);
figure; plot(soc);