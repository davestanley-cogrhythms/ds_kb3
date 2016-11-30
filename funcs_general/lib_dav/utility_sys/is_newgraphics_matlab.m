

function out = is_newgraphics_matlab

    v = version; v = str2num(v(1:3));
    out = v >= 8.4;

end