
function data_type = mode2datatype(sfc_mode)

    if sfc_mode < 20; data_type = 'SFC';
    elseif sfc_mode < 40; data_type = 'FFC';
    elseif sfc_mode >= 50 && sfc_mode <= 60; data_type = 'FR';
    else data_type = 'PSD';
    end

end
