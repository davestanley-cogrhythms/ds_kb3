
% Class containing a bunch of functions for performing swap operations on groups.
classdef group_swap
    methods (Static)

        function gr = crit(gr,swap_map_pref)
            for j = 1:length(gr)
                gr(j).criteria = swap_crit(gr(j).criteria,swap_map_pref);
            end
        end
        
        function gr = crit_move(gr,swap_map_pref)
            for j = 1:length(gr)
                gr(j).criteria = move_crit(gr(j).criteria,swap_map_pref);
            end
        end
        
        function gr = ctgs(gr,swap_map)
            for j = 1:length(gr)
                gr(j).ctgs = swap_ctgs(gr(j).ctgs,swap_map);
            end
        end
        
        function gr = cromer_shift(gr)
            for j = 1:length(gr)
                gr(j).ctgs = shift_ctgs(gr(j).ctgs); % Cheap hack to solve issue of differences in N_ctgs between cromer and roy pls, ns, ssp variables
            end
        end
        
        function gr = alt (gr)  % Adds a criteria_alt field if it does not exist. If use_alt_stage=1, then swap criteria and criteria_alt
            for j = 1:length(gr)
                    % We don't need this anymore, since all groups are now
                    % created beginning with a constructor function
%                 if ~isfield(gr(j),'criteria_alt'); gr(j).criteria_alt = []; end             % If criteria_alt field doesn't exist, create it
%                 if isempty(gr(j).criteria_alt); gr(j).criteria_alt = 2*ones(size(gr(j).criteria)); end  % If it's empty, populate it
                
%                 if use_alt_stage
                    crit_temp = gr(j).criteria_alt;
                    gr(j).criteria_alt = gr(j).criteria;
                    gr(j).criteria = crit_temp;
%                 end
                
                
            end
        end
    end
end


function crit = swap_crit(crit,swap_map)
    % Swaps criteria between origin and destination

    origin = swap_map(1,:);
    dest = swap_map(2,:);

    % swap sample and test stage criteria
    corigin = crit(:,origin);
    cdest = crit(:,dest);
    crit(:,origin) = cdest;
    crit(:,dest) = corigin; % Do dest second, incase of overwrites.
    
end


function crit = move_crit(crit,swap_map)
    % Copies criteria from origin to a temporary buffer, sets all criteria
    % to 2's, and then drops the temporary buffer into dest

    origin = swap_map(1,:);
    dest = swap_map(2,:);

    % swap sample and test stage criteria
    corigin = crit(:,origin);
    crit(:,origin) = 2;
    crit(:,dest) = corigin;

end



function ctgs = swap_ctgs (ctgs,swap_map)
%     N_ctgs_relevant = round(N_ctgs_base / 2);
%     ctgs = ctgs( ctgs < N_ctgs_relevant);

    debug_mode = 1;

    origin = swap_map(1,:);
    dest = swap_map(2,:);

    for i = 1:size(ctgs,1)      % SHould make this more efficient
        for j = 1:size(ctgs,2)
            ind = find(origin == ctgs(i,j));
            if ~isempty(ind)
                ctgs(i,j) = dest(ind);
            else
                if debug_mode warning(['Warning: No mapping found for ctgs(' num2str(i) ',' num2str(j) ') = ' num2str(ctgs(i,j))]); end
            end
        end
    end

end

function ctgs = shift_ctgs(ctgs)

    ctgs(ctgs >= 5 & ctgs <= 8) = [];
    ctgs(ctgs >= 9) = ctgs(ctgs >= 9) - 4;

end

% 
% 
% function gr = group_swap_test (gr,N_ctgs_base)
% 
%     debug_mode = 0;
% 
%     n_all = fieldnames(gr)';
%     for n = n_all
%         name = n{1};
%         if debug_mode fprintf (['Swapping to test on ' name '\n']);end
%         for j = 1:length(gr.(name))
%             gr.(name)(j).criteria = swap_crit(gr.(name)(j).criteria,N_ctgs_base);
%             gr.(name)(j).ctgs = swap_ctgs(gr.(name)(j).ctgs,N_ctgs_base);
%         end
%     end
% 
% end



