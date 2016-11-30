

function [fnum, unum] = unum2fnum2(n,fucellarray)
%     n = unum from funames1D
%     fucellarray = funames or any other cell array where each cell is a file containing vectors (representing units, electrodes, etc.)

    lengths = cellfun(@(x) length(x),fucellarray);
    clengths = cumsum(lengths);
    clengths = [0 clengths]+1; % First file starts at 1, 2nd file starts at clengths(1)+1, etc
    
    fnum=find(n >= clengths,1,'last');
    unum = n-clengths(fnum)+1;

end