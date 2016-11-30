
function Npairs = calc_Nunits_field

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
    Nelects = arrayfun(@(s) length(s.lfp_names),md);
    
    Npairs = zeros(size(Nelects));
    for k = 1:length(Nelects)
        mypairs = all_pairs(Nunits(k),Nelects(k));
        Npairs(k) = min(size(mypairs,1),Inf);  % Not sure why this Inf is here. I've dropped a warning to report if it turns up.
        if isinf(size(mypairs,1))
            warning('Inf')
        end
    end
    clear lfp1 lfp2 lfp_pairs
    
    
    
    
    %% Build string to save
    export_str = ['Nunits_field_var=('];
    for i = 1:length(Nelects)
        export_str = [export_str num2str(Npairs(i)) ' '];
    end
    export_str(end) = ')';
    
    %% Write text
    fname = 'Nunits_field.txt';
    fileID = fopen(fname,'w');
    fprintf(fileID,export_str);
    fclose(fileID);
    
    system(['chmod 775 ' fname ]);
    
end


function mypairs = all_pairs(Nunits,Nelects)
    % Generate all possible unit-electrode pairs
    [units1,lfp2] = meshgrid(1:Nunits,1:Nelects);
    mypairs = [units1(:) lfp2(:)];
    clear units1 lfp2
    Ncoherences = size(mypairs,1);
end