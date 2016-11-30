

function h = plott_matrix3D_custsolo(plotmatrix_inputs,legend_names,fb,statsmode,stats_dozscore)

    f=plotmatrix_inputs{1};
    X=plotmatrix_inputs{2};

    if statsmode == 0
        h = plott_matrix3D(plotmatrix_inputs{:});
        N_datapoints =  num2str(size(X,2));
        legend(h,[legend_names{1} ' N=' N_datapoints],[legend_names{2} ' N=' N_datapoints]);
        xlabel('Freq (Hz)');
        ylabel('SFC - normalized');
    end
        
    
    index = f >= fb(1) & f < fb(2);
    if stats_dozscore X=zscore(X); end
    dat = mean(X(index,:,:),1);
    dat = squeeze(dat);
        
    
        
        Xmu = squeeze(mean(dat,1));
        Xste = squeeze(std(dat,[],1) / sqrt(size(dat,1)));
        [h p] = ttest(diff(dat,1,2))
        
    if statsmode == 1
        
        [hbar herr ] = bar_errsig(h,Xste,Xmu);
   
        ylabel(['SFC in Range ' num2str(fb(1)) '-' num2str(fb(2)) 'Hz range']);
    elseif statsmode == 2
        nbins=20;
        subplot(211); hist(dat(:,1),nbins);
        subplot(212); hist(dat(:,2),nbins);
        h=[];
    end
       
    title(['p = ' num2str(p(1))])
    
end



