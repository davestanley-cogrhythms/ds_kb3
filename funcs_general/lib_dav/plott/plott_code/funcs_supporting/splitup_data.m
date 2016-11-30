
function Xnew = splitup_data(X,splitup)

%     FORM
%         Xnew = splitup_data(X)
%         Xnew = splitup_data(X,splitup)
%         
%     INPUTS
%         X - data matrix or vector
%         splitup - number of data points per row in Xnew
%     OUTPUTS
%         Xnew - Matrix with splitup rows and dimensionality of ndims(X)+1


    if isvector (X); X=X(:); end
    
    if nargin < 2
        splitup = round(size(X,1)/10);
    end
    
    if ndims(X) <= 2
        % Fill in the extra space at the end with zeros
        Nrows = size(X,1);
        Ncols = size(X,2);
        Nrows_new = ceil(Nrows/splitup)*splitup;
        Xnew = zeros(Nrows_new,Ncols);
        Xnew(1:Nrows,:) = X;

        Xnew = reshape(Xnew,[splitup,round(Nrows_new/splitup),Ncols]);

        clear Ntemp Ncols Nrows Nrows_new
    else
        fprintf('If using splitup, ndims of X must be <= 2 \n');
        Xnew = NaN;
        return;
    end
end