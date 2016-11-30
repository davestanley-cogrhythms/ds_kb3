
function sample_on_trials = get_sample_on_trials(sample_on_ind, ss)

    sample_on_trials = zeros(1,length(ss));
    for i = 1:length(ss)
        num_sampleons = sum(sample_on_ind(:) >= ss(i,1) & sample_on_ind(:) < ss(i,2) );
        if num_sampleons == 1
            sample_on_trials(i) = 1;
        elseif num_sampleons > 1
            sample_on_trials(i) = -1;
            fprintf('Error - sample_on_ind appears twice in the same trial!~\n');
            break;
        elseif num_sampleons < 0 
            fprintf('Should not reach this \n');
            break
        else
            % Do nothing
        end
    end

    sample_on_trials = logical(sample_on_trials);


end
