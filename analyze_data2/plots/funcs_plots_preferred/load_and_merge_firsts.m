
function [sout] = load_and_merge_firsts(file_range, sfc_mode, curr_stage_sfc,myfields)

    % Load SFC traces;
    s = loadall_ra(curr_stage_sfc,sfc_mode,myfields,[],file_range,1);
    for i = 1:length(myfields)
        myfields{i};
        sout.(myfields{i}) = s{i};
    end
    
    
%     if floor(mode) ~= 5
% 
%         
%         i=0;
%         i=i+1; sum_ctgsetli = s{i};
%         i=i+1; f =  s{i}(1,:);
%     else
%         sum_ctgsetli = [];
%         f = [];
%     end

    

end