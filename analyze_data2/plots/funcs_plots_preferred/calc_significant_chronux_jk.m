

function sig_cells = calc_significant_chronux_jk(f,Cave,Cerr1,freqval)
    
%%
    %freqval = 29;
    bh_stepup = 1;
    
    %index = find(f >= freqval, 1, 'first');
    index = find_closest(f,freqval);
    Cerr1 = Cerr1(index,:,:); Cerr1 = squeeze(Cerr1);
    Cave = Cave(index,:,:); Cave = squeeze(Cave);
    
    if ~bh_stepup
        sig_cells = Cerr1 > 0;
    else
        alpha=0.05;
        Cstd = confid_to_stdev(Cave,Cerr1,alpha);
        
        [sig_cells,pvals,pcrit] = norminv_test_bh(Cave,Cstd,alpha);
        
        
        sig_cells2=[];
        pvals2=[];
        pcrit2=[];
        for i = 1:size(Cave,2)
            [sig_cells_temp,pvals_temp,pcrit_temp] = norminv_test_bh(Cave(:,i),Cstd(:,i),alpha);
            sig_cells2 = [sig_cells2 sig_cells_temp];
            pvals2 = [pvals2 pvals_temp];
            pcrit2 = [pcrit2 pcrit_temp];
        end
        
    end
    
    sum(sig_cells(:,end));
    
end