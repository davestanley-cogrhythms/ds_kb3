

function ds = get_freq_stats(data,f,freqband_stats)

    ind = f >= freqband_stats(1) & f < freqband_stats(2);
    ds = mean(data(ind,:));
    

end