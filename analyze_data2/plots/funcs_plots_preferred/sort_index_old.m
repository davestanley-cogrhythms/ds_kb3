
function Y = sort_index_old(X,dim,I)
    % Sort data based on index returned by the sort command
    % Examples 1: Identity sort
    
%     THIS IS THE OLD VERSION OF THE CODE. LIKELY DOESN'T WORK FOR 
%     DIM = 3
    
    if dim(1) == 1
        % do nothing
    elseif dim(1) == 2
        X=X';
        I=I';
    elseif dim(1) == 3
        if length(dim) < 2 fprintf('No counting dimension specified. Assuming 2'); dim(2) = 2; end
        switch dim(2)
            case 1; X = permute(X,[3 1 2]);
            case 2; X = permute(X,[3 2 1]);
            otherwise; fprintf('Counting dimension must be 1 or 2\n'); Y=nan; return
        end
    else fprintf('Error - dim should be 1 or 2. \n');
    end
    
    [m, n] = size(X);
    for j = 1:n
        Y(:,j,:) = X(I(:,j),j,:);
    end
    
    if dim(1) == 2
        Y=Y';
    elseif dim(1) == 3
        switch dim(2)
            case 1
                Y = permute(Y,[2 3 1]);
            case 2
                Y = permute(Y,[3 2 1]);
        end
    end
    
end

