


function indhigh = find_high_beta_cells(gr,f,thresh,freqband_stats)

    
    ds = get_freq_stats(gr.data,f,freqband_stats);
    indhigh = find(ds > thresh);
    
    

    
end

