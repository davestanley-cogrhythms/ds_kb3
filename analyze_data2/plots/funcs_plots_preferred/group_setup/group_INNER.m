


function gr_out = group_INNER(group1, group2)
    % Performs inner product on group1 and group2
        % I.e. dot product (and's between element pairs of group1 and group2,
        % and then merges all of these together at the end.
    
    N1 = length(group1);
    N2 = length(group2);
    
    if N1 ~= N2
        warning('length(group1) must equal length(group2). Exiting.');
        return
    end
    
    %gr_out(1:N1*N2) = Grp;
    
    
    for i = 1:N1
        gr_out(i) = group_AND(group1(i),group2(i));
    end
    
    gr_out = group_merge(gr_out);
end
   

