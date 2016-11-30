


function sout = func_recalc_boot(buff_mode,stage_range,file_range)
   
    
    for i = 1:length(stage_range)
        sout.(stagename(stage_range(i))) = func_recalc_perms_stage(buff_mode,stage_range(i),file_range);
    end
    
    pairsdat = func_recalc_pairsmapping(buff_mode,stage_range,file_range);
    sout.pairsdat = pairsdat;
    
end

function sout = func_recalc_perms_stage(buff_mode,stage,file_range)
   
    
    %scurr = load_and_merge(file_range, buff_mode, 3, {'dCave','Cavep1','Cavep2','elapsedTime','Ntraces'});
    %scurr = load_and_merge(file_range, buff_mode, 3, {'Cave1','Cave2','phi1','phi2','Cavep1','Cavep2','phip1','phip2','elapsedTime','Ntraces'});
    scurr = load_and_merge(file_range, buff_mode, stage, {'Cave1','Cave2','phi1','phi2','ks','zerotest','dCave_bs_sig','Ntraces','Ntraces2','mypairs'});
    
    % Load "Firsts" - fields only for the first ctg
    [scurr_firsts] = load_and_merge_firsts(file_range, buff_mode, stage,{'f'});
    
    scurr_firsts.f = scurr_firsts.f(1,:);
 
    sout = catstruct(scurr, scurr_firsts);
    
end