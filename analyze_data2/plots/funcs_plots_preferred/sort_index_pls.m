


function Cavepref_out = sort_index_pls(Cave1,Cave2,I)

    test_time = 0;      % Run comparison of sorting methods to find fastest
    Cavepref = cat(4,Cave1, Cave2);    
    sz = size(Cavepref);
    I=repmat(I,[sz(1), 1, 1]);
    
    Cavepref_out=sort_index_specialized(Cavepref,4,I);
    
    
    if test_time
        tic
        Cavepref_out=sort_index(Cavepref,4,I);
        toc
        tic
        Cavepref_out=sort_index_forloops(Cavepref,4,I);     % Horribly slow, ~40 seconds
        toc
        tic
        Cavepref_out=sort_index_specialized(Cavepref,4,I);  % Seems to be fastest
        toc
    end
    
end