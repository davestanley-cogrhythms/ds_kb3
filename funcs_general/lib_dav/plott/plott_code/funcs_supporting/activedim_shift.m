

function X = activedim_shift(X,active_dim)

    Xdims = ndims(X);
    if Xdims == 3
        if active_dim == 1
            X = permute(X,[2 1 3]); % Ensures active dimension is always #2 for calculations
        elseif active_dim == 2
            X = permute(X,[1 2 3]);
        elseif active_dim == 3
            X = permute(X,[1 3 2]);
            %X = permute(X,[1 active_dim 5-active_dim]); % Swaps dim2 and dim3 if active_dim = 3
        else
            fprintf('Error. Active dim must equal 1, 2, or 3. exiting\n');
            return;
        end
    end
    if Xdims == 2
        if active_dim == 1 || active_dim == 2
            X = permute(X,[3-active_dim active_dim]);   % Ensures active dimension is always #2 for calculations
        elseif active_dim == 3
            X = permute(X,[1 3 2]);
        else
            fprintf('Error. Active dim must equal 1 or 2. exiting\n');
            return;
        end


end