

function s = loadall_units(return_field,file_list,file_range)
    % loadall_units(stage,mode,return_field,file_list,file_range)
    % all inputs optional
    % file_list - list of files to load (default, loads all in path)
    % file_range - range of file numbers within the list to load (probably can never use this)

    if ~exist('return_field','var'); return_field = 'unit_ind'; end
    %if ~exist('return_field','var'); return_field = 'unit_array'; end

    datapath = getpath('path_metadata');
    
    if ~exist('file_list','var')
        file_list = get_filelist(datapath);
        
    end
    
    if ~exist('file_range','var')
        file_range = 1:length(file_list);
        file_range = 1:3;
    end
    
    md = cellfun(@load,fullfile(datapath,file_list(file_range)));
    unitsca=horzcat({md(:).(return_field)});
    
    files_examined = file_list(file_range);
    units_examined=horzcat({md(:).unit_names});
    units_examined=convert_unit_underscores(units_examined);
    units_per_file = cellfun(@length,units_examined);
    
    s=v2struct(unitsca,files_examined,units_examined,units_per_file);
    
end



