

function Y = sort_index_specialized(X,dim,I)
    % Sort data based on index returned by the sort command
    % CAn handle matrices
    
    % For the following description, let dim = 1 and [i,j,k] = size(X)
    % Performs exactly the same operation as sort_index, except it is
    % optimized for the case when j and k are large and the sorting
    % dimension, i, is small.
    % It's likely going to be slower when i is larger.
    % Examples: See below

    if dim(1) == 1
        % do nothing
    elseif dim(1) == 2
        X=permute(X,[2 1 3]);
        I=permute(I,[2 1 3]);
    elseif dim(1) == 3
        X = permute(X,[3 1 2]);
        I = permute(I,[3 1 2]);    
    elseif dim(1) == 4
        X = permute(X,[4 1 2 3]);
        I = permute(I,[4 1 2 3]);    
    else fprintf('Error - dim should be between 1 and 3. \n');
    end
    
    % If I contains some singleton dimensions, expand.
    sz = size(X);
    sz2 = size(I);
    sze = sz - sz2; % Size by which to expand X
    sze = sze + 1;  % This way, dimensions that are the same size will be 1; dimensions that are larger will be rep-matted as needed
    I = repmat(I,sze);

    
    
    index_vals = sort(unique(I(:))); % Find all index values
    Y = zeros(size(X));
    for i = 1:length(index_vals);
        Y(i,:) = X(I == i);
    end
    
        
    if dim(1) == 2
        Y=permute(Y,[2 1 3]);
    elseif dim(1) == 3
        Y = permute(Y,[2 3 1]);
    elseif dim(1) == 4
        Y = permute(Y,[2 3 4 1]);
    end
    
end

% 
% 
% % % % Examples for testing sort_index.m - These are actually the tests I used to debug the code
% clear
% % % Generate 3D Matrix
% A=1:4; A=A';
% B=1:3;
% C=1:2; C = permute(C,[1 3 2]);
% 
% B=B*10; C=C*100;
% A=repmat(A,[1,3,2]);
% B=repmat(B,[4,1,2]);
% C=repmat(C,[4,3,1]);
% 
% X_orig=A+B+C;
% X=X_orig;
% 
% I1a=[4 3 2 1]; I1a=I1a(:);
% sort_index(X,1,I1a)
% 
% I1b = [3 2 1];
% sort_index(X,2,I1b)
% 
% I1c = [2 1]; I1c = permute(I1c,[1 3 2]);
% sort_index(X,3,I1c)
% 
% 
% I2a = [4 3 2 1; 1 2 3 4; 1 2 3 4]; I2a = I2a'
% sort_index(X,1,I2a)
% 
% I2b = [4 3 2 1; 1 2 3 4]'; I2b = permute(I2b,[1,3,2])
% sort_index(X,1,I2b)
% 
% % Tests the same conditions that are used for running SFC
% I2c = [2 1 1; 1 2 2]; I2c = permute(I2c,[3 2 1])
% sort_index(X,3,I2c)
% 
% I3 = [1 2 3 4]'; I3 = repmat(I3,[1,3,2]); I3(:,1,1) = [4:-1:1]
% sort_index(X,1,I3)
% 
% 
% 
% 
% % % Generate 2D Matrix
% X=X_orig(:,:,1); X=X-100
% 
% % Reverse along dimension one
% I1a=[4 3 2 1]; I1a=I1a(:)
% sort_index(X,1,I1a)
% 
% % Reverse along dimension two
% I1b = [3 2 1]
% sort_index(X,2,I1b)
% 
% %Tests the same conditions that are used for sorting neuron firing rates (src)
% I1b = [3 2 1; 1 2 3; 1 2 3; 1 2 3]
% sort_index(X,2,I1b)
% 
% 
% 
% 
% 
% 
