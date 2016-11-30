

function [bad_indices, badfiles] = list_badfiles(mystat, threshold, lfp_names)
    

    % Use forloops for readability and laziness, since this code is fast
    % anyawys
    bad_indices = {};
    badfiles = {};
    k=1;
    
    for i = 1:length(mystat)
        for j = 1:length(mystat{i})
            
            
%             % This code I used it as a short hack to allow a less strict threshold level for Monkey L than for Monkey O
%             if i < 41
%                 threshold = 5; else
%                 threshold = get_clipping_threshold;
%             end
            
            bad_indices{i}(j) = mystat{i}(j) > threshold;
            if bad_indices{i}(j)
                badfiles{k} = lfp_names{i}{j};
                k=k+1;
            end
        end
    end
end

