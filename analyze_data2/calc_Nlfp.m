
function Nunits = calc_Nlfp

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
    Nunits = arrayfun(@(s) min(length(s.lfp_names),Inf),md);  % Not sure why this Inf is here. I've dropped a warning to report if it turns up.
    Nunits2 = arrayfun(@(s) length(s.lfp_names),md);
    if any(isinf(Nunits2)); warning('Inf detected'); end; 
    
    
    %% Build string to save
    export_str = ['Nlfp_var=('];
    for i = 1:length(Nunits)
        export_str = [export_str num2str(Nunits(i)) ' '];
    end
    export_str(end) = ')';
    
    %% Write text
    fname = 'Nlfp.txt';
    fileID = fopen(fname,'w');
    fprintf(fileID,export_str);
    fclose(fileID);
    
    system(['chmod 775 ' fname ]);
    
end
