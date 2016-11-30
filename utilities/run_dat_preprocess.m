


extract_metadata_and_ints
extract_metadata_merge
extract_lfp_sample_delay
identify_clipping(-2);
identify_clipping(-1);
identify_clipping(0);
identify_clipping(1);
identify_clipping(2);
identify_clipping(3);
identify_clipping(4);
identify_clipping(5);
identify_clipping(6);
identify_clipping(7);
identify_clipping(8);
identify_60Hz

filter60_lfp_sample_delay % Not implemented for ad6 (ad6 was using unfiltered data)

count_percent_switch_trials;
morphs_classify;    % Classifies the images in the morphs folder (i.e. catdog psycho testing) as cat/dog/fat/thin using SVM.
                    % Then converts condition codes in conds files to
                    % say which morphs are used for test images.
morphs_organize;    % This takes all the morph classifications (which are specific
                    % to the conds files (of which there are extra, 1:84) and
                    % maps them to the file names of interest (i.e. the
                    % metadata file names - 1:79)