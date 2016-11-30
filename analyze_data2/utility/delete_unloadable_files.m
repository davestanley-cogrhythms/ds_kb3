

function delete_unloadable_files(sfc_mode, curr_stage, folder_prefix)
    % Finds and deletes any files of size 0 bytes in i0j0_split folder


    %% Setup vars
    if ~exist('sfc_mode'); sfc_mode = 22.4018111; end
    if ~exist('curr_stage'); curr_stage = 2; end
    if ~exist('folder_prefix'); folder_prefix = 'test2'; end
    
    %% First, get list of unmerged files
    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix, ~, ~, ~, ~, ~, permutation_test] = build_sfcmode(sfc_mode, mode_subgroups);
    outpath_base = fullfile(getpath('path_buffer_curr'),[folder_prefix 'mode_' num2str(sfc_mode,12) fname_suffix]);
    stage_base = ['stage_' num2str(curr_stage)];
    outpath = fullfile(outpath_base,stage_base,'i0j0_split');
    
    %% Code to list all files

    % Suspected zeros
    d1 = dir(fullfile(outpath,'*.mat'));
    n1 = {d1.name};
    sz = [d1.bytes];

    %% Estimate the unloadable files based on file size
    clear missing_files 
    j=1;
    missing = sz < (median(sz)-1*std(sz)); % All small files - 1 standard deviation below median
    missing_ind = find(missing);
    missing_files = n1(missing);
    
    
    %% Confirm detected files can't be loaded
    missing2 = false(size(missing));
    for i = 1:length(missing_ind)
        try
            load(fullfile(outpath,missing_files{i}));
        catch
            missing2(missing_ind(i)) = true;
        end
    end
    
    
    %% Update missing file list
    missing_ind = find(missing2);
    missing_files = n1(missing2);
    
    
    %% Write zeros to file
    save(fullfile(outpath_base,['unloadable_ind_' stage_base '.txt']),'missing_ind','-ascii');
    fileID = fopen(fullfile(outpath_base,['unloadable_files_' stage_base '.txt']),'w');
    for i = 1:length(missing_ind)
        fprintf(fileID,[missing_files{i} '\n']);
    end
    fclose(fileID);

    %% Move to deletion
    mkdir(fullfile(outpath,'deleteme'));
    for i = 1:length(missing_ind)
    %     i
        %system(['rm ' missing_files{i}]);
        system(['mv ' fullfile(outpath,missing_files{i}) ' ' fullfile(outpath,'deleteme')],'-echo');
    end

    %% Make sure all files empty
    d2 = dir(fullfile(outpath,'deleteme'));
    figure; plot([d2.bytes]); ylabel('Size (bytes)'); xlabel('Files in deleteme');

%     %% Delete folder
%     system(['rm -rf ' fullfile(outpath,'deleteme')]);

end