

function logiout = build_logical(array, N)
    %% Builds an array of logicals from an array of indices
    
%     tic
    out = false(1,N);
    out(array) = true;
    logiout = out;
%     toc 
    
   
%     % Old code
%     tic
%     out = logical(zeros(1,N));
%     for i = array
%         out(i) = 1;
%     end
%     
%     logiout2 = out;
%     toc
end