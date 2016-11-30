
function [pairstype, pairs_mode] = mode2pairstype(sfc_mode)

    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    pairs_mode = mode_subgroups(1);
    
    switch pairs_mode
        case 0
            pairstype='u0';
        case 2
            pairstype='u0adj';
        case 3
            pairstype='ue';
        case 4
            pairstype='ee';
        case 5
            pairstype='uu';
        case 6
            pairstype='e0';
        case 7
            pairstype='u0';
    end
    
end
