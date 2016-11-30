


function gr_out = group_OUTER(group1, group2)
    % Performs outer product on group1 and group2
    
    N1 = length(group1);
    N2 = length(group2);
    
    %gr_out(1:N1*N2) = Grp;
    
    k=0;
    for i = 1:N1
        for j = 1:N2
            k=k+1;
            gr_out(k) = group_AND(group1(i),group2(j));
        end
    end
    

end


