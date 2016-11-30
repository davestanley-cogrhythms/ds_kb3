
function currspike_ctg_out = thin_spks(currspike_ctg, prob_thinning)

    sz = size(currspike_ctg);

    ind = find(currspike_ctg);
    ind = ind(randperm(length(ind)));
    
    tot_spks = length(ind);
    num_spks = tot_spks*prob_thinning;
    chosen_spks = ind(1:num_spks);
    
    currspike_ctg_out = zeros(sz);
    currspike_ctg_out(chosen_spks) = 1;

end

