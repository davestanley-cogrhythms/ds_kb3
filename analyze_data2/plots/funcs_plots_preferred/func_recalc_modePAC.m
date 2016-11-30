
function adat = func_recalc_modePAC (spect_mode,stage_range,file_range)

    for i = 1:length(stage_range)

        extract_field = {'Cave','p'};
            
        % Load fields for all ctgs
        [scurr] = load_and_merge(file_range, spect_mode, stage_range(i), extract_field);
        
        
        % Load "Firsts" - fields only for the first ctg
        %[scurr_firsts] = load_and_merge_firsts(file_range, spect_mode, 5, {'T','F','sum_ctgsetli','N_ctgs_base','N_ctgs_extras','Nwind'});      % params field still not present in some data
        myfields_first = {'sum_ctgsetli','f','f2','N_ctgs_base','N_ctgs_extras','params'};
        [scurr_firsts] = load_and_merge_firsts(file_range, spect_mode, stage_range(i), myfields_first);
        
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
        sz = size(scurr.Cave);
        scurr.Cave = reshape(scurr.Cave,[Nrows,Ncols,sz(2),sz(3)]);
        scurr.p = reshape(scurr.p,[Nrows,Ncols,sz(2),sz(3)]);
        
        % Do the permute
        scurr.Cave = permute(scurr.Cave,[1,3,4,2]);
        scurr.p= permute(scurr.p,[1,3,4,2]);
        
        % Merge
        scurr = catstruct(scurr, scurr_firsts);
        
        adat.(stagename(stage_range(i))) = scurr;
        
    end
end


function y = take_middles(x)
    %y = repmat(x,[2,1]);
    y = [x(1:end-1); x(2:end)];
    y = mean(y,1);
end