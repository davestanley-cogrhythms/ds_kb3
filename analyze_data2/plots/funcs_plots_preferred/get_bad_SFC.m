


function b = get_bad_SFC(X)
    
    X = X(:,:,:);
    Xn = isnan(X) | isinf(X);
    Xn = sum(Xn,1);
    Xn = sum(Xn,3);
    Xn = squeeze(Xn);
    
    b = Xn >= 1;
    
    
   

end