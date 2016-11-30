
function [sens_merged,funames1D] = loadall_morph_sens(stage)

    morph_mode = 2;

    % Load morph sensitivities
    morphpath = fullfile(getpath('path_buffer_curr'),'morphs');
    filename = ['mode_' num2str(morph_mode) '_stage_' num2str(stage)];
    morph_name = fullfile(morphpath,filename);
    ms = load(morph_name);
    
    sens_merged = cat(1,ms.sens{:});
    funames1D = cat(2,ms.funames{:});
    
end

