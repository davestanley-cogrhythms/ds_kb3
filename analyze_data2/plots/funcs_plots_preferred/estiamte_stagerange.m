

function stage_range = estiamte_stagerange(sfc_mode)

    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix] = build_sfcmode(sfc_mode, mode_subgroups);
    
    datapath = fullfile(getpath('path_buffer_curr'),['mode_' num2str(sfc_mode,12) fname_suffix]);
    
    out = dir(fullfile(datapath,'stage*'));
    foldernames = {out.name};
    
    
    stage_range = zeros(length(foldernames),1);
    for i = 1:length(foldernames)
        curr_folder = foldernames{i};
        stage_range(i) = str2num(curr_folder(7:end));
    end
    
    if isempty(stage_range)
        error(['No stages found. Path searched: ' datapath]);
    end


end