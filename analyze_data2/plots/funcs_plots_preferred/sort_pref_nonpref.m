

function [pref,nonpref,Cavepref] = sort_pref_nonpref(f,pls1,pls2,freqband_stats)

    plot_debug = 0;
        data_type = [];
        var_str = ['Coherence'];
    
    ind = find_closest(f,mean(freqband_stats));
    Cavepref = cat(3,pls1(:,:), pls2(:,:));
    Cavef = Cavepref(ind,:,:);
    [~, I]=sort([Cavef],3,'descend');
    sz = size(Cavepref);
    I=repmat(I,[sz(1), 1, 1]);
    Cavepref=sort_index(Cavepref,3,I);
    if plot_debug
        figure; plott_matrix3D(Cavepref(:,:,:));xlim([0 100]); xlabel('Freq (Hz)'); ylabel(var_str); legend(['Preferred ' data_type],['Non-preferred ' data_type]);
        title(['Preferred vs Non-Preferred']);
    end

    pref = Cavepref(:,:,1);
    nonpref = Cavepref(:,:,2);


end