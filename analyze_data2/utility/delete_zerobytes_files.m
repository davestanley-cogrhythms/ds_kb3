

function delete_zerobytes_files(sfc_mode, curr_stage, folder_prefix)
    % Finds and deletes any files of size 0 bytes in i0j0_split folder


    %% Setup vars
    if ~exist('sfc_mode'); sfc_mode = 22.4718111; end
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

    %% Search
    clear missing_files 
    j=1;
    missing_ind = find(sz == 0);
    for i = 1:length(missing_ind);
        missing_files{i} = n1{missing_ind(i)};
    end
    
    %% Write zeros to file
    save(fullfile(outpath_base,['zerobytes_ind_' stage_base '.txt']),'missing_ind','-ascii');
    fileID = fopen(fullfile(outpath_base,['zerobytes_files_' stage_base '.txt']),'w');
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

 
    %% Delete folder
    system(['rm -rf ' fullfile(outpath,'deleteme')]);

end