

function Z = diff_pad(X,n,dim,padval)

    if nargin < 4; padval=0;end
    if nargin < 3; dim=1; end
    if nargin < 2; n=[]; end
    
    if isvector(X); X=X(:); end
    
    Y = diff(X,n,dim);
    
    % Get the number of rows to add
    Norig = size(X,1);
    Nfin = size(Y,1);
    sizdiff = Norig - Nfin;
    
    % Split these among the top and bottom
    num_rows = floor(sizdiff/2); padstart = add_padding(X,num_rows,padval);
    num_rows = ceil(sizdiff/2); padfin = add_padding(X,num_rows,padval);
    
    Z = cat(1,padstart,Y,padfin);
    
end

function padding = add_padding(X,num_rows,padval)
    matdims=size(X);
    matdims(1) = num_rows;
    padding = padval * ones(matdims);

end

