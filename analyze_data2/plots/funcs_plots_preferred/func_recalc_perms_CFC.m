


function sout = func_recalc_perms(buff_mode,stage_range,file_range)
   
    
    for i = 1:length(stage_range)
        sout.(stagename(stage_range(i))) = func_recalc_perms_stage(buff_mode,stage_range(i),file_range);
    end
    
end

function sout = func_recalc_perms_stage(buff_mode,stage,file_range)
   
    
    %scurr = load_and_merge(file_range, buff_mode, 3, {'dCave','Cavep1','Cavep2','elapsedTime','Ntraces'});
    %scurr = load_and_merge(file_range, buff_mode, 3, {'Cave1','Cave2','phi1','phi2','Cavep1','Cavep2','phip1','phip2','elapsedTime','Ntraces'});
    scurr = load_and_merge(file_range, buff_mode, stage, {'Cave1','Cave2','pvals_Cave','Ntraces','Ntraces2','mypairs'});
    
    % Load "Firsts" - fields only for the first ctg
    [scurr_firsts] = load_and_merge_firsts(file_range, buff_mode, 3,{'f','f2'});
    
    scurr_firsts.f = scurr_firsts.f(1,:);
    scurr_firsts.f2 = scurr_firsts.f2(1,:);
    scurr_firsts.forig = scurr_firsts.f;
    scurr_firsts.f2orig = scurr_firsts.f2;
    scurr_firsts.f = take_middles(scurr_firsts.f);
    scurr_firsts.f2 = take_middles(scurr_firsts.f2);
    
    % Now unpack Cave by reshaping it into a matrix. Then permute it to
    % be in the format that the rest of the code expects
    Nrows = length(scurr_firsts.f);
    Ncols = length(scurr_firsts.f2);
    sz = size(scurr.Cave1);
    scurr.Cave1 = reshape(scurr.Cave1,[Nrows,Ncols,sz(2)]);
    scurr.Cave2 = reshape(scurr.Cave2,[Nrows,Ncols,sz(2)]);
    scurr.pvals_Cave = reshape(scurr.pvals_Cave,[Nrows,Ncols,sz(2)]);

    % Do the permute
    scurr.Cave1 = permute(scurr.Cave1,[1,3,4,2]);
    scurr.Cave2 = permute(scurr.Cave2,[1,3,4,2]);
    scurr.pvals_Cave = permute(scurr.pvals_Cave,[1,3,4,2]);

 
    sout = catstruct(scurr, scurr_firsts);
    
end



function y = take_middles(x)
    %y = repmat(x,[2,1]);
    y = [x(1:end-1); x(2:end)];
    y = mean(y,1);
end