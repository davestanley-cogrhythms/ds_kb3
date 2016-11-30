
function [currlfp_ctg_out, currspike_ctg_out] = thin_spks_bootstrap(currlfp_ctg, currspike_ctg, prob_thinning, Nbootstraps)

    retaining_factor = 0.8;  % The final total number of spikes in currlfp_ctg_out should be approximately
                             % retaining_factor * the original number of
                             % spikes. When this is 1.0, we try to use all
                             % the data, although there may be some
                             % redundancy

    Ntrials = size(currspike_ctg,2);
    sz = size(currspike_ctg);
    Ntrials_max = 10000;
    
    % Estimate # bootstraps
    if ~exist('Nbootstraps','var')
        [Nbootstraps,Ntrials_out] = est_Nbootstraps(prob_thinning, sz(2), retaining_factor, Ntrials_max);

%         Nbootstraps_old = floor(1/prob_thinning * (retaining_factor));   % If we keep only 10% of spikes on each bootstrap, we want to conduct 10 bootstrapping trials to bring
%                                                                      % us back to 100.
%         % Prune back # bootstraps if # trials gets too high (if this is too high, things get really slow)
%         Ntrials_out = min(sz(2)*Nbootstraps_old,Ntrials_max);
%         Nbootstraps_old = floor(Ntrials_out / sz(2));
%         
%         Ntrials_out_old = min(sz(2)*Nbootstraps_old);
    else
        Ntrials_out = min(sz(2)*Nbootstraps);
    end
    
    
    currlfp_ctg_out = zeros(sz(1),Ntrials_out);
    currspike_ctg_out = currlfp_ctg_out;
    
    for i = 1:Nbootstraps

        ind = find(currspike_ctg);
        ind = ind(randperm(length(ind)));

        tot_spks = length(ind);
        num_spks = tot_spks*prob_thinning;
        chosen_spks = ind(1:min(floor(num_spks),length(ind)));
        %chosen_spks = ind(1:floor(num_spks));

        currspike_ctg_temp = zeros(sz);
        currspike_ctg_temp(chosen_spks) = 1;
        
        currspike_ctg_out(:,(i-1)*Ntrials+1:i*Ntrials) = currspike_ctg_temp;
        
        if ~isempty(currlfp_ctg)
            currlfp_ctg_out(:,(i-1)*Ntrials+1:i*Ntrials) = currlfp_ctg;
        end
    
    end

end

