
function adat = func_recalc_sfc(stage_range,file_range,sfc_mode)

    if ~exist('stage_range','var')
       stage_range = [];
    end
    if isempty(stage_range)
        stage_range = estiamte_stagerange(sfc_mode);
    end
 
    % Load SFC for the current stage in stage_range.
    for i = 1:length(stage_range)
        % Load fields for all ctgs
        
        myfields = {'Cave','phiave','phistd','Ntraces_below_thresh','Ntraces','spikerate_mu','p','confC','sum_ctgsetli_good','sum_ctgsetli_maxes','mypairs'};
        
        [scurr] = load_and_merge(file_range, sfc_mode, stage_range(i),myfields);
        
        %if ~isempty(scurr.spikerate_mu)        % Might want to use this clause if ever combine this code with func_recalc_pairs
            scurr.spikerate_mu = scurr.spikerate_mu(:,end); scurr.spikerate_mu = scurr.spikerate_mu(:)';
        %end
        
        
        % Load "Firsts" - fields only for the first ctg
        myfields_first = {'sum_ctgsetli','f','f2','N_ctgs_base','N_ctgs_extras','params'};
        [scurr_firsts] = load_and_merge_firsts(file_range, sfc_mode, stage_range(i),myfields_first);
        scurr_firsts.f = scurr_firsts.f(1,:);
        if ~isempty(scurr_firsts.f2); scurr_firsts.f2 = scurr_firsts.f2(1,:); end
        scurr_firsts.N_ctgs_base = scurr_firsts.N_ctgs_base(1);
        scurr_firsts.N_ctgs_extras = scurr_firsts.N_ctgs_extras(1);
        
        % Remove N_curr* fields - we will inherit these from query functions
        % instead
%         scurr_firsts = rmfield(scurr_firsts, 'N_ctgs_base');
%         scurr_firsts = rmfield(scurr_firsts, 'N_ctgs_extras');
        
        % Merge
        scurr = catstruct(scurr, scurr_firsts);
        
        % Add to struct
        adat.(stagename(stage_range(i))) = scurr;
    end

    pairsdat = func_recalc_pairsmapping(sfc_mode,stage_range,file_range);
    adat.pairsdat = pairsdat;
    
end