

function s = loadall_ra(stage,sfc_mode,return_field,file_list,file_range,first_cell_only)
    % Loads all data output by run_analysis.m
    % data_organized has the format {ctg,filenum}(1,unitnum,f)

    if ~exist('stage','var'); stage = 3; end
    if ~exist('sfc_mode','var'); sfc_mode = 2; end
    if ~exist('return_field','var'); return_field = 'Cave'; end
    if ~exist('first_cell_only','var'); first_cell_only = 0; end
    
    %if ~exist('return_field','var'); return_field = 'S1ave'; end
    %if ~exist('return_field','var'); return_field = 'S2ave'; end


    % Check if we're in .2 mode - "Adjacent electrode" mode or switch
    % trials mode. Adjust folder suffix name appropriately.
    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix] = build_sfcmode(sfc_mode, mode_subgroups);
    
    datapath = fullfile(getpath('path_buffer_curr'),['mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(stage)]);
    
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
    
    if isempty(file_range)
        file_range = 1:length(file_list);
        %file_range = 1:3;
    end
    
    data_default = cellfun(@load,fullfilec(datapath,file_list(file_range)));
    if ~first_cell_only
        if iscell(return_field)
            for i = 1:length(return_field)
                s{i} = format_data(data_default, return_field{i});
            end
        else
            s = format_data(data_default, return_field);
        end
    else
        
        if iscell(return_field)
            for i = 1:length(return_field)
                s{i} = format_data_firstcell(data_default, return_field{i});
            end
        else
            s = format_data_firstcell(data_default, return_field);
        end
    end

end



function s = format_data(data_default, return_field)

    datacell = vertcat(data_default(:).sfc);          % Pull out sfc field
    datacell = permute(datacell,[2 1]);       % Take transpose - this will make it easier to stack subsequent files with cell2mat
    
    % Squeeze out any cells not containing return field
    emptycells = cellfun(@(s) ~isfield(s,return_field),datacell);
    datacellsq = datacell(all(~emptycells,2),:);

    %if isfield(datacellsq{1},return_field);
        b = cellfun(@(s) s.(return_field),datacellsq,'UniformOutput',0); % Dynamic fieldnames
    %else
    %    b = cellfun(@(s) [],datacellsq,'UniformOutput',0);
    %end
    data_organized = cellfun(@(s) permute(s,[ndims(s)+1, 2, 1, 3:ndims(s)]),b,'UniformOutput',0); % Align the data. This ensures the 1st dimension
                                                                                                  % is empty, the 2nd dimension is cells, 3rd...ndims dimensions
                                                                                                  % are extra data (times, frequencies, etc)
    
    files_examined = cellfun(@(s) s.('filename'),datacell(1,:),'UniformOutput',0);
    units_examined = cellfun(@(s) s.('unit_names'),datacell(1,:),'UniformOutput',0);
    units_examined=convert_unit_underscores(units_examined);
    file_unit_names = fname_uname_merge (files_examined,units_examined);
    units_per_file = cellfun(@length,units_examined);
    
    data = data_organized;
    s=v2struct(data,files_examined,units_examined,units_per_file,file_unit_names);
    
end


function s = format_data_firstcell(data_default, return_field)
    % This function takes out a field of the first entry in sfc structure
    % cell array (sfc{1}). This field must be a row matrix. It replicates
    % it corresponding to the number of neurons in that file and returns
    % the result.
    % 


    % Get number of cells per file
    Ncells = arrayfun(@(s) length(s.sfc{end}.Ntraces),data_default);
    
    datacell = vertcat(data_default(:).sfc);          % Pull out sfc field
    datacell = datacell(:,1);                         % Take only the first column, corresponding to sfc{1}
    
    if isfield(datacell{1},return_field);
        b = cellfunu(@(s) s.(return_field),datacell);     % Extract the chosen field
        b = cellfunu(@(v) force_vector_horizontal(v), b);
    else
        b = cellfunu(@(s) [],datacell);     % Field doesn't exist - return empty
    end
    
    
    bmat = cellfunu(@(x,N) repmat(x,N,1),b,num2cell(Ncells(:)));
    
    s = cell2mat(bmat);
end


function v = force_vector_horizontal(v)
    if isvector(v)
        v = v(:)';
    end
end

