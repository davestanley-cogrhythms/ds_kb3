

function md = loadall_md(file_list,file_range)
    % loadall_units(stage,mode,return_field,file_list,file_range)
    % all inputs optional
    % file_list - list of files to load (default, loads all in path)
    % file_range - range of file numbers within the list to load (probably can never use this)

    datapath = getpath('path_metadata');
    
    if ~exist('file_list','var')
        file_list = get_filelist(datapath);
    end
    
    if isempty(file_list)
        file_list = get_filelist(datapath);
    end
    
    if ~exist('file_range','var')
        file_range = 1:length(file_list);
        %file_range = 1:3;
    end
    
    md = cellfun(@load,fullfilec(datapath,file_list(file_range)));
    %md = cellfun(@load, fullfilec( );
    
end




