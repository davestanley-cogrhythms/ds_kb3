

function [bad_clip, bad_60,fnames_clip] = loadall_bad_data(curr_stage_badclipping,file_range,md)

    curr_stage_badchannels_60Hz = 5;
    unit_LFP_pairs = get_allunit_LFP_mapping(md);

    % Load bad data & merge
    path_badchannels = fullfile(getpath('path_badchannels'));
    name_clip = ['id_clipping_s' num2str(curr_stage_badclipping) '.mat'];
    name_60Hz = ['id_60Hz_s' num2str(curr_stage_badchannels_60Hz) '.mat'];
    s_clip = load(fullfile(path_badchannels,name_clip));
    s_60Hz = load(fullfile(path_badchannels,name_60Hz));
    li_clip = cat(2,s_clip.clippingstat{:});
    li_60 = cat(2,s_60Hz.linestat{:});
    % li_hist = cat(2,s_clip.Xint_hist{:});
    % li_psd = cat(2,s_60Hz.psdout{:});
    % li_f = s_60Hz.f;
    bad_clip_li = cat(2,s_clip.bad_indices{:});
    bad_60_li = cat(2,s_60Hz.bad_indices{:});

    % Optional - use custom thresholds for selecting bad files
    clth = get_clipping_threshold_trials;
    bad_clip_li = li_clip > clth(2);
    % bad_60_li = li_60 > get_60Hz_threshold;
    % bad_clip_li = li_clip > 1.5;
    % bad_60_li = li_60 > 0.1;
    
    bad_clip = bad_clip_li(unit_LFP_pairs(2,:));
    bad_60 = bad_60_li(unit_LFP_pairs(2,:));

    % % Import filenames (optional)
    fnames_clip = s_clip.fnames;
    flnames_clip = s_clip.flnames;
    
    flnames_clip1D = cat(2,flnames_clip{:});
    
    fnames_60 = s_60Hz.fnames;
    flnames_60 = s_60Hz.flnames;
    flnames_601D = cat(2,flnames_60{:});
    
    if sum(~strcmp(fnames_clip,fnames_60)) ~= 0; fprintf(['File mismatch! Likely some files are missing for either '...
        'numspikes or sfc. Aborting \n']); end

    if sum(~strcmp(flnames_clip1D,flnames_601D)) ~= 0; fprintf(['File mismatch! Likely some files are missing for either '...
        'numspikes or sfc. Aborting \n']); end

    [a b]  = unum2fnum(100,fnames_clip,flnames_clip);          % Restore the original file and lfp numbers from the merged list, for testing
    

end