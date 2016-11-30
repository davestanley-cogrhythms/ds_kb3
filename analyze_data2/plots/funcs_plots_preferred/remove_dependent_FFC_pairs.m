

function goodout = remove_dependent_FFC_pairs(mypairs,chosen0)

    %%
    if nargin < 2
        chosen0 = true(size(mypairs,1),1);
    end
    
%     tic
    
    selection_mode = 1; % 1-draw random electrode; 2-draw random pair
                        % These have similar processing times.
    
    % Resize chosen
    chosen = chosen0(:);
    
    % This is our pool of unique electrodes that are available
    mypairs2 = mypairs(chosen,:);
    ue = unique(mypairs2(:));
    
    % Initialize goodout - list of mypairs that, together, contain
    % only unique electrodes.
    goodout = false(size(mypairs,1),1);
    while ~isempty(ue)
        
        if selection_mode == 1
            % Draw a random electrode.
            [ue_curr] = draw_random(ue);

            % Find all pairs that contain this electrode. Pairs must be
            % eligible (i.e. they cannot be amongst chosen)
            limatch1 = (mypairs(:,1) == ue_curr | mypairs(:,2) == ue_curr) & chosen(:);

            % Pick one pair at random to keep for goodout
            pcurr = draw_random(find(limatch1));
            goodout(pcurr) = true;

            % Find the complementary electrode
            pair = mypairs(pcurr,:);
            if pair(1) == ue_curr; ue_comp = pair(2);
            elseif pair(2) == ue_curr; ue_comp = pair(1);
            else
                warning('should not reach this'); keyboard;
            end

            % Find all pairs that contain complement
            limatch2 = (mypairs(:,1) == ue_comp| mypairs(:,2) == ue_comp) & chosen(:);

        elseif selection_mode == 2
            
            % Draw a random pair from set of eligibles
            [pcurr] = draw_random(find(chosen));
            
            % Mark as good
            goodout(pcurr) = true;
            
            % Get current and complementary electrodes
            ue_curr = mypairs(pcurr,1);
            ue_comp = mypairs(pcurr,2);
            
            limatch1 = (mypairs(:,1) == ue_curr | mypairs(:,2) == ue_curr) & chosen(:);
            limatch2 = (mypairs(:,1) == ue_comp| mypairs(:,2) == ue_comp) & chosen(:);
            
            
        end
        
        % Remove both from chosen
        chosen = chosen & ~(limatch1 | limatch2);
        
        
        % Reassess list of chosen electrodes
        mypairs2 = mypairs(chosen,:);
        ue = unique(mypairs2(:));
        
    end
    
    

%     toc
end


function [y] = draw_random(x)
    % Draws a random element y from the set x

    y = x(unidrnd(length(x),1,1));

end

