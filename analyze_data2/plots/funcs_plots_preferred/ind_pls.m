
function I = ind_pls(Cave1,Cave2,f,freqband_stats)

    ind = find_closest(f,mean(freqband_stats));
    Cavepref = cat(4,Cave1, Cave2);             % Use dimension 4 just incase it's a 3D matrix.
    Cavef = Cavepref(ind,:,:,:);
    [~, I]=sort([Cavef],4,'descend');
end

