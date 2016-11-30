
function Nunits = calc_Nunits

    %%Setup paths
    currdir = pwd;
    cd ..
    setup_paths
    setup_paths
    cd (currdir);
    
    %% load md
    path_md = fullfile(getpath('path_metadata'),'sansunits');
    file_list = get_filelist(path_md);
    md = cellfun(@(f) load (fullfile(path_md,f)),file_list);     % Load metadata
    Nunits = arrayfun(@(s) min(length(s.unit_names),Inf),md);  % Not sure why this Inf is here. I've dropped a warning to report if it turns up.
    Nunits2 = arrayfun(@(s) length(s.unit_names),md);
    if any(isinf(Nunits2)); warning('Inf detected'); end; 
    
    
    %% Build string to save
    export_str = ['Nunits_var=('];
    for i = 1:length(Nunits)
        export_str = [export_str num2str(Nunits(i)) ' '];
    end
    export_str(end) = ')';
    
    %% Write text
    fname = 'Nunits.txt';
    fileID = fopen(fname,'w');
    fprintf(fileID,export_str);
    fclose(fileID);
    
    system(['chmod 775 ' fname ]);
    
end
