
function out = and_across_rows(in)
    Nrows = size(in,2);
    out=sum(in,2) == Nrows ;

end
