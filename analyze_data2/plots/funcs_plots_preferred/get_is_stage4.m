

function out = get_is_stage4(sfc_mode)
    
    sfc_mode_group = floor(sfc_mode);
    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    
    if get_is_spectrogram(sfc_mode) || ...              % If it's a spectrogram, then we're in mode 4
            (sfc_mode_group == 52 && mode_subgroups(4) > 0)  % If it's a windowed time series, then we're in mode 4
        out = 1;
    else
        out = 0;
    end
    

end