
function [file_list] = get_filelist(data_path,wildcard)

    if nargin < 2
        wildcard = '*.mat';
    end

    %Nfiles = 79;
    randomize_files = 0;


    %data_path = getpath('path_roydata_orig');
    % fullfile(filesep,'media','davestanley','Dave','Users','Computer','data','CRC_Miller_data');
    % other useful commands - filesep, fileparts, pathsep
    % /media/davestanley/Dave/Users/Computer/data/CRC_Miller_data

    available_files = dir(fullfile(data_path,wildcard));
    if randomize_files
        select_files = randperm(length(available_files));
    else
        select_files = 1:length(available_files);
    end
    
    clear file_list;
    Nfiles = length(available_files);       % Use maximum number of available files
    for i = 1:Nfiles
        file_list{i} = available_files(select_files(i)).name;
    end


end