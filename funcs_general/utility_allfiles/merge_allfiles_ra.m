function [ra_dat, ra_fname, ra_funame] = merge_allfiles_ra(ra_struct)
    % Pulls the data out of the structure returned by loadall_ra. This data
    % is organized based on file names. This script converts it into per
    % units. This function handles conditions, and therefore should only
    % be used with the ouput of loadall_ra / run_analysis
    
    % INPUT
    % Input is the output structure of loadall_ra.m
    % OUTPUT
    % Output data is a matrix of the form (x,y,z), where x is the data
    % vector, y is the unit, and z is the category scheme (1-9).
    % Also outputs file names and file-unit names cell arrays

    plot_on = 0;
    
    % Pull out data from structure
    ra_dat=cell2mat(ra_struct.data);
    ra_dat = permute(ra_dat,[3,2,1, 4:ndims(ra_dat)]);
    
    % Pull out filenames
    ra_fname=ra_struct(:).files_examined;
    
    
    % Pull out unit names
    ra_funame=ra_struct(:).file_unit_names; ra_funame=horzcat(ra_funame{:});
    
    if plot_on
        figure; plott_ani(ra_dat);
    end
    
end


