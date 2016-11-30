function pls_sort = get_pls_sort(out_pls, out_sort, sort_on, do_sort_irr, freqband_stats)
        
    Cave1 = out_sort.Cave1;
    Cave2 = out_sort.Cave2;

    % Sort morphs
    if sort_on
        [IA] = ind_pls(Cave1(:,:,1),Cave2(:,:,1),out_sort.abscissa,freqband_stats);  % Sch A
        [IB] = ind_pls(Cave1(:,:,2),Cave2(:,:,2),out_sort.abscissa,freqband_stats);  % Sch B
        if do_sort_irr
            [IAirr] = ind_pls(Cave1(:,:,3),Cave2(:,:,3),out_sort.abscissa,freqband_stats);  % Sch A
            [IBirr] = ind_pls(Cave1(:,:,4),Cave2(:,:,4),out_sort.abscissa,freqband_stats);  % Sch B
        else
            IAirr = IA;
            IBirr = IB;
            warning('Maybe redo for irrelevant deciders.');
        end
    else
        IA = [];
        IB = [];
        IAirr = [];
        IBirr = [];
    end
    

    sz = size(out_pls.pls);
    pls_sort = zeros([sz(1),sz(2),28]);

    del_p=8;
    del_ps = 7;
    % Sch A Relevant [1,2]; arrange into "Roy" type structure.
    pls_sort(:,:,[3+0*del_ps,5+0*del_ps]) = sort_index_pls_wrapper(out_pls.pls(:,:,1+0*del_p),out_pls.pls(:,:,2+0*del_p),IA,sort_on);  % 60%
    pls_sort(:,:,[2+0*del_ps,6+0*del_ps]) = sort_index_pls_wrapper(out_pls.pls(:,:,1+1*del_p),out_pls.pls(:,:,2+1*del_p),IA,sort_on);  % 80%
    pls_sort(:,:,[1+0*del_ps,7+0*del_ps]) = sort_index_pls_wrapper(out_pls.pls(:,:,1+2*del_p),out_pls.pls(:,:,2+2*del_p),IA,sort_on);  % 100%
    pls_sort(:,:,[4+0*del_ps]) = out_pls.pls(:,:,del_p*3+1);

    % Sch B Relevant [3,4]
    pls_sort(:,:,[3+1*del_ps,5+1*del_ps]) = sort_index_pls_wrapper(out_pls.pls(:,:,3+0*del_p),out_pls.pls(:,:,4+0*del_p),IB,sort_on);  % 60%
    pls_sort(:,:,[2+1*del_ps,6+1*del_ps]) = sort_index_pls_wrapper(out_pls.pls(:,:,3+1*del_p),out_pls.pls(:,:,4+1*del_p),IB,sort_on);  % 80%
    pls_sort(:,:,[1+1*del_ps,7+1*del_ps]) = sort_index_pls_wrapper(out_pls.pls(:,:,3+2*del_p),out_pls.pls(:,:,4+2*del_p),IB,sort_on);  % 100%
    pls_sort(:,:,[4+1*del_ps]) = out_pls.pls(:,:,del_p*3+2);


    % Sch A IrrRelevant (i.e. cued for doing fat/thin, but pls is measuring
    % cat/dog
    pls_sort(:,:,[3+2*del_ps,5+2*del_ps]) = sort_index_pls_wrapper(out_pls.pls(:,:,5+0*del_p),out_pls.pls(:,:,6+0*del_p),IAirr,sort_on);  % 60%
    pls_sort(:,:,[2+2*del_ps,6+2*del_ps]) = sort_index_pls_wrapper(out_pls.pls(:,:,5+1*del_p),out_pls.pls(:,:,6+1*del_p),IAirr,sort_on);  % 80%
    pls_sort(:,:,[1+2*del_ps,7+2*del_ps]) = sort_index_pls_wrapper(out_pls.pls(:,:,5+2*del_p),out_pls.pls(:,:,6+2*del_p),IAirr,sort_on);  % 100%
    pls_sort(:,:,[4+2*del_ps]) = out_pls.pls(:,:,del_p*3+4);        % Sch B cue + 50% irrelevant

    % Sch B IrrRelevant (i.e. cued for doing cat/dog, but pls is measuring
    % fat/thin
    pls_sort(:,:,[3+3*del_ps,5+3*del_ps]) = sort_index_pls_wrapper(out_pls.pls(:,:,7+0*del_p),out_pls.pls(:,:,8+0*del_p),IBirr,sort_on);  % 60%
    pls_sort(:,:,[2+3*del_ps,6+3*del_ps]) = sort_index_pls_wrapper(out_pls.pls(:,:,7+1*del_p),out_pls.pls(:,:,8+1*del_p),IBirr,sort_on);  % 80%
    pls_sort(:,:,[1+3*del_ps,7+3*del_ps]) = sort_index_pls_wrapper(out_pls.pls(:,:,7+2*del_p),out_pls.pls(:,:,8+2*del_p),IBirr,sort_on);  % 100%
    pls_sort(:,:,[4+3*del_ps]) = out_pls.pls(:,:,del_p*3+3);        % Sch A cue + 50% irrelevant
    
end


function Cavepref_out = sort_index_pls_wrapper(Cave1,Cave2,I,sort_on)
    
    if sort_on
        Cavepref_out = sort_index_pls(Cave1,Cave2,I);
    else
        
        Cavepref_out = cat(4,Cave1, Cave2);

    end

end

