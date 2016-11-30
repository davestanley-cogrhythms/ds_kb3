
function [curr_fname, curr_fieldname] = get_matfiles_names(path_matfiles_store, prefix, sfc_mode)

    curr_fieldname = [prefix mode2modename(sfc_mode)];
    curr_fname = fullfile(path_matfiles_store,[curr_fieldname '.mat']);
    
end