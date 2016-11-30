
function [bad_switch] = get_bad_switch(sum_ctgsetli, mode_subgroups, isolate_clustered_protocols)


Ncells = size(sum_ctgsetli,1);
bad_switch = false(1,Ncells);                   % Take all cells

if mode_subgroups(2) == 1                       % This means we ARE using switch trials
    N_ctgs = size(sum_ctgsetli,2) - 2;                  % There are N-2 categories. N-1 is all switch trials; N is all trials
    sw_ind = N_ctgs + 1;
    all_ind = N_ctgs + 2;
    if isolate_clustered_protocols == 1         % Take only cells with clustering of switch trials
        bad_switch = (sum_ctgsetli(:,sw_ind) ./ sum_ctgsetli(:,all_ind) * 100) > 25;
        bad_switch = bad_switch(:)';
    elseif isolate_clustered_protocols == 2     % Take only cells with no clustering of switch trials
        bad_switch = (sum_ctgsetli(:,sw_ind) ./ sum_ctgsetli(:,all_ind) * 100) <= 25;
        bad_switch = bad_switch(:)';
    end
end


end