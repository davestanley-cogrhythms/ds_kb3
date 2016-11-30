
function [Nbootstraps,Ntrials_out] = est_Nbootstraps(prob_thinning, Ntrials_orig, retaining_factor, Ntrials_max)


    if nargin < 4
        Ntrials_max = 10000;
    end
    
    if nargin < 3
        retaining_factor = 0.8;  % The final total number of spikes in currlfp_ctg_out should be approximately
                             % retaining_factor * the original number of
                             % spikes. When this is 1.0, we try to use all
                             % the data, although there may be some
                             % redundancy
    end
    

    
    


    Nbootstraps = floor(1/prob_thinning * (retaining_factor));   % If we keep only 10% of spikes on each bootstrap, we want to conduct 10 bootstrapping trials to bring
                                                                 % us back to 100.
	Nbootstraps = max(1,Nbootstraps);                            % Cannot be less than 1!

    % Prune back # bootstraps if # trials gets too high (if this is too high, things get really slow)
    Ntrials_out = min(Ntrials_orig*Nbootstraps,Ntrials_max);        % Can't have more than Ntrials_max
    Nbootstraps = floor(Ntrials_out / Ntrials_orig);                % Recalculate appropriate Nbootstraps based on this
    Ntrials_out = min(Ntrials_orig*Nbootstraps);                    % Re-estimate Ntrials_out


    
end

