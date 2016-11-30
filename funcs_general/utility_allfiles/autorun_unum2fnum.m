


function [fnum, unum, curr_funame] = autorun_unum2fnum(n)
    % Returns file and unit numbers corresponding to n, without requiring
    % you to actually load any data. This just makes some assumptions and
    % handles loading the data itself.

    curr_stage = 3;
    mode = 2.2;
    
    file_range = 1:79;
    sfc_struct = loadall_ra(curr_stage,mode,{'Cave'},[],file_range);
    
    fnames = sfc_struct{1}.files_examined;
    funames = sfc_struct{1}.file_unit_names;
    
    [fnum, unum, curr_funame] = unum2fnum(n,fnames,funames);

end