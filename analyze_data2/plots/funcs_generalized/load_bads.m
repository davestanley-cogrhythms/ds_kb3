function [bad_any,bad_clip,bad_60,fnames_clip] = load_bads(sfc_mode,curr_stage_badclipping,opts_exclude,md,mypairs,funames1D)

    % Load file range
    [~, file_range] = estimate_filelist(sfc_mode, curr_stage_badclipping);

    % Pull out ue_pairs
    [~, sfc_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix, do_adjacent, ctgsetli_mode, thinning_mode, tapers_mode, baseline_subtract, permutation_test, ue_pairs, coh_debias_mode, do_partial_coherence, ctgsetli_mode2] = build_sfcmode(sfc_mode, sfc_subgroups);

    % Pull in excludes
    exclude_clipping = opts_exclude.exclude_clipping;
    exclude_60 = opts_exclude.exclude_60;
    exclude_nans = opts_exclude.exclude_nans;
    excludeL = opts_exclude.excludeL; 
    excludeO = opts_exclude.excludeO; 
    remove_dependent = opts_exclude.remove_dependent;       % Remove dependent electrode pairs for FFC
   

    % Load bad files
    if any(ue_pairs == [0,2])           % SFC
        [bad_clip, bad_60,fnames_clip] = loadall_bad_data(curr_stage_badclipping,file_range,md);
    elseif any(ue_pairs == [4,6,3])     % FFC, PSD, or SFC with units vs elects
        [bad_clip, bad_60,fnames_clip] = loadall_bad_data_lfp(curr_stage_badclipping,file_range,md);
    elseif ue_pairs == 7;
        [bad_clip, bad_60,fnames_clip] = loadall_bad_data(curr_stage_badclipping,file_range,md);
        bad_clip = false(size(bad_clip));
        bad_60 = false(size(bad_60)); % All good!
    end


    % Get bad files
    bad_any = false(size(bad_clip));
    if exclude_clipping; bad_any = bad_any | bad_clip; end
    if exclude_60; bad_any = bad_any | bad_60; end

    % If in FFC mode, convert bads from single to paired
    if ue_pairs == 4
        bad_any = any([bad_any(mypairs(:,1))',bad_any(mypairs(:,2))'],2)';
        bad_60 = any([bad_60(mypairs(:,1))',bad_60(mypairs(:,2))'],2)';
        bad_clip = any([bad_clip(mypairs(:,1))',bad_clip(mypairs(:,2))'],2)';
    elseif ue_pairs == 3
        bad_any = bad_any(mypairs(:,2));
        bad_60 = bad_60(mypairs(:,2));
        bad_clip = bad_clip(mypairs(:,2));
        % Mark as bad all pairs that have units on the matching electrode
        warning('Need to remove units on same electrode!');
    end

    ismonkeyL_temp = get_ismonkeyL(funames1D);
    if excludeL; bad_any = bad_any | ismonkeyL_temp; end
    if excludeO; bad_any = bad_any | (~ismonkeyL_temp); end
    if remove_dependent; bad_any = ~remove_dependent_FFC_pairs(mypairs,~bad_any)'; end
    clear ismonkeyL_temp


end