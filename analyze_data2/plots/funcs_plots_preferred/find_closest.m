


function indc = find_closest(x,val)

    % Finds the index in vector x that is closest to value val. This can be
    % either larger or smaller.
    % 
    % This is a more advanced way of just doing:
    %    index = find(x >= val,1,'first');
    % which will only find closest index >= val
    do_comparison=0;
    
    if isnan(val) || isempty(val); error('val is not properly assigned'); end
    
    xd = abs(x-val);
    [~,indc] = min(xd);

    if do_comparison
        index = find(x >= val,1,'first');
        index2 = find(x <= val,1,'last');
    end

end