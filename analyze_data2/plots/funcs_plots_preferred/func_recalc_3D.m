
function scurr = func_recalc_3D (spect_mode,file_range)

        extract_field = {'Ssp','S2AVE','ConfC'};
            
        % Load fields for all ctgs
        [scurr] = load_and_merge(file_range, spect_mode, 5, extract_field);
        
        % Load "Firsts" - fields only for the first ctg
        %[scurr_firsts] = load_and_merge_firsts(file_range, spect_mode, 5, {'T','F','sum_ctgsetli','N_ctgs_base','N_ctgs_extras','Nwind'});      % params field still not present in some data
        [scurr_firsts] = load_and_merge_firsts(file_range, spect_mode, 5, {'T','F','sum_ctgsetli','N_ctgs_base','N_ctgs_extras','params','Nwind','F2'});
        
        % Merge
        scurr = catstruct(scurr, scurr_firsts);
        
end

