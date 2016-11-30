

function s = func_recalc_md(file_range)



buff_filelist = get_filelist(fullfile(getpath('path_metadata'),'sansunits'));

if nargin < 1
    md = cellfun(@(filename) load (fullfile(getpath('path_metadata'),'sansunits',filename)),buff_filelist(:)');
else
    md = cellfun(@(filename) load (fullfile(getpath('path_metadata'),'sansunits',filename)),buff_filelist(file_range));
end


s.md = md;

end