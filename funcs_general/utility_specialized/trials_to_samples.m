

function samplesli = trials_to_samples( chosen_trialsli, sample_on_trials)

%     % Old code
%     chosen_trialsli = chosen_trialsli(:);
%     sample_on_trials = sample_on_trials(:);
%     
%     Ntrials = length(chosen_trialsli);
%     Nsamples = sum(sample_on_trials);
%     
%     ns = (1:Nsamples)';
%     
%     trial_sample_mapper = zeros(Ntrials,1);
%     trial_sample_mapper(logical(sample_on_trials)) = ns;
%     
%     samples = trial_sample_mapper(logical(chosen_trialsli));
%     
%     samplesli = build_logical(samples,Nsamples);
    
    
    % Faster way
    chosen_trialsli = chosen_trialsli(:);
    sample_on_trials = sample_on_trials(:);
    
    samplesli = chosen_trialsli(sample_on_trials);
    samplesli=samplesli(:)';                        % Make sure output is a row - same format as before
    
    
end

