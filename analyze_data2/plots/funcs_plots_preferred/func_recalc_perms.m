


function adat = func_recalc_perms(sfc_mode,stage_range,file_range)

    if ~exist('stage_range','var')
       stage_range = [];
    end
    if isempty(stage_range)
        stage_range = estiamte_stagerange(sfc_mode);
    end

    for i = 1:length(stage_range)
        adat.(stagename(stage_range(i))) = func_recalc_perms_stage(sfc_mode,stage_range(i),file_range);
    end

    pairsdat = func_recalc_pairsmapping(sfc_mode,stage_range,file_range);
    adat.pairsdat = pairsdat;
    
end

function adat = func_recalc_perms_stage(sfc_mode,stage,file_range)
   
    
    %scurr = load_and_merge(file_range, sfc_mode, 3, {'dCave','Cavep1','Cavep2','elapsedTime','Ntraces'});
    %scurr = load_and_merge(file_range, sfc_mode, 3, {'Cave1','Cave2','phi1','phi2','Cavep1','Cavep2','phip1','phip2','elapsedTime','Ntraces'});
    scurr = load_and_merge(file_range, sfc_mode, stage, {'Cave1','Cave2','phi1','phi2','pvals_Cave','pvals_phi','Ntraces','Ntraces2','mypairs','zsc_Cave','zsc_phi'});
    %scurr = load_and_merge(file_range, sfc_mode, stage,
    %{'Cave1','Cave2','phi1','phi2','pvals_Cave','pvals_phi','Ntraces','Ntraces2','mypairs','sum_ctgsetli_good1','sum_ctgsetli_maxes1','sum_ctgsetli_good2','sum_ctgsetli_maxes2','zsc_Cave','zsc_phi'}); %  This one is the most recent, but good1/good2 are redundant with Ntraces/Ntraces2
    
    % Load "Firsts" - fields only for the first ctg
    [scurr_firsts] = load_and_merge_firsts(file_range, sfc_mode, stage,{'f','f2','sum_ctgsetli','stagesir'});
    scurr_firsts.f = scurr_firsts.f(1,:);
    if ~isempty(scurr_firsts.f2); scurr_firsts.f2 = scurr_firsts.f2(1,:); end
    scurr_firsts.stagesir = scurr_firsts.stagesir(1,:);
    
 
    adat = catstruct(scurr, scurr_firsts);
    
end