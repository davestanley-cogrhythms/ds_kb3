

function sig_cells = calc_significant_chronux(f,C,confC,freqval)
    

    %index = find(f >= freqval, 1, 'first');
    index = find_closest(f,freqval);
    
    C = C(index,:,:);
    C = squeeze(C);
    sig_cells = C > confC;
    
end