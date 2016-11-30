
function [x2,y2] = reduce_xy_plusminus(x,y,split_plusminus)
    switch split_plusminus
        case 2
            ind = x >= 0;   % Only take positives
        case 3
            ind = x < 0;    % Only take negatives
        otherwise
            ind = true(size(x));    % Take all
    end
    
    x2 = x(ind);
    y2 = y(ind);
end