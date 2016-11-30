



function lfp_all = merge_allfiles_lfp(lfp)
% function lfp_all = merge_allfiles_lfp(lfp)
% This function merges the cell array "lfp" containing all files, each with multiple channels
% into one long cell array. Should generally be used after the following command
% lfp = cellfun(@(filename) load (fullfile(path_lfp_sample_delay,filename)), file_list(file_range)); % Load all files. This could take a while...


    lfp_all = permute(lfp,[1 3 2]);
    lfp_all = cat(3, {lfp_all(:).lfp_sample});
    lfp_all = cellfunu(@electrodes_to_cell_array,lfp_all);
    lfp_all = cat(2,lfp_all{:});
    
    
end
    

function cell_out = electrodes_to_cell_array (mat_in)
    sz = size(mat_in);
    cell_out = mat2cell(mat_in,sz(1),sz(2),ones(1,sz(3)));
    cell_out = cell_out(:)';
end