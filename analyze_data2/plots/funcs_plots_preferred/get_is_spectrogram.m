

function out = get_is_spectrogram(sfc_mode)
    sfc_mode_group = floor(sfc_mode);
    if any(sfc_mode_group == [3, 23, 33, 45])
        out=1;
    else
        out=0;
    end

end