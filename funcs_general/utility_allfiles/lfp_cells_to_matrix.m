
function lfp_mat = lfp_cells_to_matrix(lfp_all)

    lfp_mat = cellfunu(@(x) x(:), lfp_all);
    lfp_lengths = cellfun(@length, lfp_mat);
    minlen = min(lfp_lengths);
    lfp_mat = cellfunu(@(x) x(1:minlen),lfp_mat);
    lfp_mat = cell2mat(lfp_mat);
    
end
    