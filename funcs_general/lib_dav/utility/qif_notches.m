

function out= qif_notches(in,width)
    in = in(:);
    x = in-width/2;
    y = in+width/2;
    out = [x y];
end

