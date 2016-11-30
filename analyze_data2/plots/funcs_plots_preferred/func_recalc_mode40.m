
function adat = func_recalc_mode40 (spect_mode,stage_range,file_range)

    for i = 1:length(stage_range)

        extract_field = {'Cave','Ptr','Cerr1','Cerr2'};
            
        % Load fields for all ctgs
        [scurr] = load_and_merge(file_range, spect_mode, stage_range(i), extract_field);
        
        scurr.Ptr = permute(scurr.Ptr,[1,4,2,3]);
        scurr.Cave = permute(scurr.Cave,[1,4,2,3]);
        scurr.Cerr1 = permute(scurr.Cerr1,[1,4,2,3]);
        scurr.Cerr2 = permute(scurr.Cerr2,[1,4,2,3]);
        
        % Load "Firsts" - fields only for the first ctg
        %[scurr_firsts] = load_and_merge_firsts(file_range, spect_mode, 5, {'T','F','sum_ctgsetli','N_ctgs_base','N_ctgs_extras','Nwind'});      % params field still not present in some data
        myfields_first = {'sum_ctgsetli','f','N_ctgs_base','N_ctgs_extras','params'};
        [scurr_firsts] = load_and_merge_firsts(file_range, spect_mode, stage_range(i), myfields_first);
        
        % Merge
        scurr = catstruct(scurr, scurr_firsts);
        
        adat.(stagename(stage_range(i))) = scurr;
        
    end
end


