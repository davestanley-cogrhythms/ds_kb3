function make_sample_on_per_trial
  %% save out to path defined in my_paths via query defined in funcs_general/funcs_queries
  outpath = fullfile(getpath('path_sample_on_times_per_trial'));
  
  %% get all session filenames
  file_list = dir(fullfile(getpath('path_lfp'),'*.mat'));
  
  for i = 1:length(file_list)
    disp(['Processing file ' num2str(i) ' of ' num2str(length(file_list)) ': ' file_list(i).name]);
    make_sample_on_times_per_trial(file_list(i).name,outpath);
  end

end

function make_sample_on_times_per_trial(filename,outpath)  
  load (fullfile(getpath('path_metadata'),filename)); % Load metadata
  %% need sample_on_times and sample_on_trials

  num_trials = length(sample_on_trials);
  sample_on_times_per_trial = zeros(num_trials, 1, 'double');

  n = 1;
  for i = 1:num_trials
    %% if we have a sample onset time, put it in our new vector
    if sample_on_trials(i)
      sample_on_times_per_trial(i) = sample_on_times(n);
      n = n + 1;
    else
      %% if not, put NAN
      sample_on_times_per_trial(i) = NaN;
    end
  end
save(fullfile(outpath,filename),'sample_on_times_per_trial');

end