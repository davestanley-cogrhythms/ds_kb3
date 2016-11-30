

function extract_metadata_and_ints(extract_metadata,extract_intmatrix)
    % Extracts metadata and integer lfp matrices from Jefferson's original
    % mat files. Saves them to the default paths (specified by getpath_roy_data)
    % extract_metadata = 1 or 0 to extract metadata
    % extract_intmatrix = 1 or 0 to extract and save integer matrices
    clc
    if ~exist('extract_metadata','var'); extract_metadata = 1; end
    if ~exist('extract_intmatrix','var'); extract_intmatrix = 1; end
    
    % Operates in several modes:
    % extract_metadata - extracts and saves metadata from Jefferson's original files
    % extract_intmatrix - extracts and saves data as integers from Jefferson's original files
    
    
    data_path = getpath('path_roydata_orig');
    metadata_outputpath = fullfile(data_path,'metadata');
    intmatrix_outputpath = fullfile(data_path,'intmatrix');
    
    if extract_metadata
        mkdir_dave (metadata_outputpath);
    end
    if extract_intmatrix
        mkdir_dave (intmatrix_outputpath);
    end
    
    file_list = dir([data_path filesep '*.mat']);
    
    parfor i = 1:length(file_list)
            % This is really messy, but basically just passes the whole
            % workspace to this function. Used for parallel processing.
        extract_singlefile(data_path,file_list,extract_metadata,extract_intmatrix,i,metadata_outputpath,intmatrix_outputpath)
    end
    
end


function extract_singlefile(data_path,file_list,extract_metadata,extract_intmatrix,i,metadata_outputpath,intmatrix_outputpath)

    fprintf ('Processing file %d of %d, name %s \n',i,length(file_list),file_list(i).name);
    s = load (fullfile(data_path,file_list(i).name));

    if extract_metadata
        s2 = rmfield(s,'lfp_array');
        save (fullfile(metadata_outputpath,file_list(i).name),'-struct','s2');
    end

    if extract_intmatrix
        lfpmat_int16 = int16(cell2mat(s.lfp_array));
        save (fullfile(intmatrix_outputpath,file_list(i).name),'lfpmat_int16');
    end

        

end