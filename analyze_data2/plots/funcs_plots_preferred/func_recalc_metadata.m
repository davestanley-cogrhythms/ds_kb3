    
function func_recalc_metadata(filename,file_range)
    % Load Metadata
    file_list = get_filelist(fullfile(getpath('path_metadata'),'sansunits'));
    md = cellfun(@(filename) load (fullfile(getpath('path_metadata'),'sansunits',filename)),file_list(file_range));
    
    % Load evoked stuff
    ev_struct = loadall_ra(5,9,'E',[],file_range);
    [ev, ev_fnames, ev_funames1D] = merge_allfiles_ra(ev_struct); clear ev_struct
    evdat = v2struct(ev_fnames, ev_funames1D,ev);
    
    
    save(filename,'md','evdat');
end