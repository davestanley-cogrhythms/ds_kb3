
function [s] = load_and_merge(file_range, sfc_mode, curr_stage_sfc,myfields)

    % Load SFC traces
    sfc_struct = loadall_ra(curr_stage_sfc,sfc_mode,myfields,[],file_range);
    i=0;
    for i = 1:length(myfields)
        [s.(myfields{i}), s.fnames, s.funames1D] = merge_allfiles_ra(sfc_struct{i});
        s.(myfields{i}) = squeeze(s.(myfields{i}));
    end
    

end