function y = fisher_trans_shifted(x)
    
    % Identity (range of fisher transform is 0 to inf)
    a=1;
    b=0;
    
    % To stretch from -1 to 1 before doing fisher transform (range is -inf
    % to inf)
    a=2;
    b=-1;
    
    y = atanh(a*x+b);

end