function x = fisher_trans_inv(y)

    a=1;
    b=0;
    %y = atanh(a*x-b);

    x = (tanh(y)+b) / a;
    
end