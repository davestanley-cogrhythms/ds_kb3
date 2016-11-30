
function Npairs = calc_Nunits_units

    %%Setup paths
    currdir = pwd;
    cd ..
    setup_paths
    setup_paths
    cd (currdir);
    addpath(genpath('./funcs_supporting_local'));
    
    %% load md
    path_md = fullfile(getpath('path_metadata'),'sansunits');
    file_list = get_filelist(path_md);
    md = cellfun(@(f) load (fullfile(path_md,f)),file_list);     % Load metadata
    Nunits = arrayfun(@(s) length(s.unit_names),md);
    
    
    Npairs = zeros(size(Nunits));
    for k = 1:length(Nunits)
        mypairs = upper_triangle_pairs(Nunits(k));
        Npairs(k) = min(size(mypairs,1),Inf);  % Not sure why this Inf is here. I've dropped a warning to report if it turns up.
        if isinf(size(mypairs,1))
            warning('Inf');
        end
    end
    clear lfp1 lfp2 lfp_pairs
    
    
    
    
    %% Build string to save
    export_str = ['Nunits_units_var=('];
    for i = 1:length(Nunits)
        export_str = [export_str num2str(Npairs(i)) ' '];
    end
    export_str(end) = ')';
    
    %% Write text
    fname = 'Nunits_units.txt';
    fileID = fopen(fname,'w');
    fprintf(fileID,export_str);
    fclose(fileID);
    
    system(['chmod 775 ' fname ]);
    
end
