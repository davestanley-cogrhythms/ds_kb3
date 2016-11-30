

function count_i0j0_split_files(sfc_mode, curr_stage, folder_prefix)
    % Creates a text document that tells number of files in
    % each i0j0_split folder
    % Also lists the files that are missing, if any


    %% Setup vars
    if ~exist('sfc_mode'); sfc_mode = 22.4718111; end
    if ~exist('curr_stage'); curr_stage = 2; end
    if ~exist('folder_prefix'); folder_prefix = 'test2'; end
    
    %% First, get list of i0j0_split files
    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix, ~, ~, ~, ~, ~, permutation_test] = build_sfcmode(sfc_mode, mode_subgroups);
    outpath_base = fullfile(getpath('path_buffer_curr'),[folder_prefix 'mode_' num2str(sfc_mode,12) fname_suffix]);
    stage_base = ['stage_' num2str(curr_stage)];
    outpath = fullfile(outpath_base,stage_base,'i0j0_split');
    
    %% Code to list all files

    % Suspected zeros
    d1 = dir(fullfile(outpath,'*.mat'));
    n1 = {d1.name};
    
    %% Search for missing files
    load('example_full_list_lfp_pairs.mat');    % Example full filename list of lfp pairs for Roy (loads n2 and d2)
                                                % Based on find_missing_files.m
    clear missing_files 
    missing_files = {};
    missing_ind = [];
    j=1;
    for i = 1:length(n2)
        temp=strcmp(n1,n2{i});
        if sum(temp) == 0;
            missing_files{j} = n2{i};
            missing_ind = [missing_ind i];
            j=j+1;
        end

    end


    %% Write data to textfile
    fileID = fopen(fullfile(outpath_base,['numfiles_' stage_base '.txt']),'w'); % Open for writing, discard existing contents, if any
    
    % Write number of files
    temp = length(n1);
    fprintf(fileID,[num2str(temp) '\n']);
    
    % Write missing indices
    fprintf(fileID,'-- Missing indices -- \n');
    for i = 1:length(missing_ind)
        fprintf(fileID,[num2str(missing_ind(i)) '\n']);
    end
    
    % Write missing filenames
    fprintf(fileID,'-- Missing files -- \n');
    for i = 1:length(missing_files)
        fprintf(fileID,[missing_files{i} '\n']);
    end
   
    
    fclose(fileID);

end