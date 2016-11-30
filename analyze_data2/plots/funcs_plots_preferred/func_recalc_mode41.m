
function adat = func_recalc_mode41 (sfc_mode,stage_range,file_range)

    if ~exist('stage_range','var')
       stage_range = [];
    end
    if isempty(stage_range)
        stage_range = estiamte_stagerange(sfc_mode);
    end

    for i = 1:length(stage_range)

        extract_field = {'Cave','Cerr1','Cerr2','mypairs'};
            
        % Load fields for all ctgs
        [scurr] = load_and_merge(file_range, sfc_mode, stage_range(i), extract_field);
        

        % Load "Firsts" - fields only for the first ctg
        %[scurr_firsts] = load_and_merge_firsts(file_range, sfc_mode, 5, {'T','F','sum_ctgsetli','N_ctgs_base','N_ctgs_extras','Nwind'});      % params field still not present in some data
        myfields_first = {'sum_ctgsetli','f','f2','params','sum_ctgsetli','stagesir'};
        [scurr_firsts] = load_and_merge_firsts(file_range, sfc_mode, stage_range(i), myfields_first);
        scurr_firsts.f = scurr_firsts.f(1,:);
        if ~isempty(scurr_firsts.f2); scurr_firsts.f2 = scurr_firsts.f2(1,:); end
        
        % Merge
        scurr = catstruct(scurr, scurr_firsts);
        
        adat.(stagename(stage_range(i))) = scurr;
        
    end
    
    pairsdat = func_recalc_pairsmapping(sfc_mode,stage_range,file_range);
    adat.pairsdat = pairsdat;
end


