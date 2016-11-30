

function convert_preferred_units



path_metadata = fullfile(getpath('path_metadata'));

file_list = get_filelist(path_metadata);
Nfiles = length(file_list);


tot_Nunits=  0;
for i = 1:Nfiles
    fprintf ('Counting units in file %d of %d, name %s \n',i,length(file_list),file_list{i});
    filename = file_list{i};
    md = load (fullfile(getpath('path_metadata'),filename));     % Load metadata
    [mdout{i} md_matrix{i}] = calc_metadata_preferred(filename,md);
end



end


