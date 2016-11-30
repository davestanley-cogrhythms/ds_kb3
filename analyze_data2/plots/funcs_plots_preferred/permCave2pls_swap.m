function [pls_new] = permCave2pls_swap(s_perm,perm2pls_dophi)
    %% permCave2pls_swap

    f1='Cave1';
    f2='Cave2';
    if perm2pls_dophi
        f1='phi1';
        f2='phi2';
    end

    % Query size of 1st 4 dimensions
        % This guarantees to return entries even if size(x,3) = 1
    sz1 = arrayfun(@(x) size(s_perm.(f1),x),1:4);
    sz2 = arrayfun(@(x) size(s_perm.(f2),x),1:4);
    
    N = sz1(3) + sz2(3);
    sz_new = sz1; sz_new(3)=N;
    pls_new = zeros(sz_new);
    
    pls_new(:,:,1:2:N) = s_perm.(f1);   % Do Cave1 in odd positions
    pls_new(:,:,2:2:N) = s_perm.(f2);   % Do Cave2

end