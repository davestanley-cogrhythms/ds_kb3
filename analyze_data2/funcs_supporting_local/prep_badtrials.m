
function bad_trials = prep_badtrials(curr_stage_badchannels)

    
    % Load bad data & merge
    path_badchannels = fullfile(getpath('path_badchannels'));
    name_clip = ['id_clipping_s' num2str(curr_stage_badchannels) '.mat'];
    
    s_clip = load(fullfile(path_badchannels,name_clip));
    bad_trials = s_clip.bad_trials;

end
