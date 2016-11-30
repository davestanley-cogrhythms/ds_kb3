

function [fname_suffix, do_adjacent, ctgsetli_mode, thinning_mode, tapers_mode, baseline_subtract, permutation_test, ue_pairs, coh_debias_mode, do_partial_coherence, ctgsetli_mode2, Nwind_mode] = build_sfcmode(sfc_mode, mode_subgroups)

    sfc_mode_subgroup = mode_subgroups(1);
    sfc_mode_subgroup2 = mode_subgroups(2);
    sfc_mode_subgroup3 = mode_subgroups(3);
    sfc_mode_subgroup4 =  mode_subgroups(4);
    sfc_mode_subgroup5 =  mode_subgroups(5);
    sfc_mode_subgroup6 =  mode_subgroups(6);

    if sfc_mode_subgroup == 2; do_adjacent = 1; % Look at adjacent electrode
    else do_adjacent = 0; end
    
    ue_pairs = sfc_mode_subgroup;
    
    ctgsetli_mode = sfc_mode_subgroup2;    % Include in ctgsetli 
    
    thinning_mode = sfc_mode_subgroup3;
    
    tapers_mode = sfc_mode_subgroup4;               % 0-Default tapers; 1-long tapers; 2-long tapers + Jackknife
                                                    % 3-Tapers=[1,1]; 4-Tapers=[3,5]; 5-Tapers=[5,9]
    baseline_subtract = sfc_mode_subgroup5;
    
    permutation_test = sfc_mode_subgroup6;
    
    coh_debias_mode = mode_subgroups(7);
    
    do_partial_coherence = mode_subgroups(8);
    
    ctgsetli_mode2 = mode_subgroups(9); % Second digit in ctgsetli_mode
    
    Nwind_mode = mode_subgroups(10);    % Size of window for spectrogram
    
    switch sfc_mode_subgroup
        case 0; adj_mode = '';
        case 2; adj_mode = '_adj';
        case 3; adj_mode = '_ue';
        case 4; adj_mode = '_ee';
        case 5; adj_mode = '_uu';
        case 7; adj_mode = '_u0';   % Unit only
        otherwise; adj_mode = '';
    end
    
    
    
    
    if floor(sfc_mode) ~= 5
        switch ctgsetli_mode2
            case {0,2}      % Default 
                cmt = ctgsetli_mode_2_string(ctgsetli_mode,ctgsetli_mode2);
            case {1,3}      % For merging relevant and irrelevant
                cmt = ctgsetli_mode_2_string(ctgsetli_mode,ctgsetli_mode2);
                cmt = [cmt 'Mrg'];
            otherwise
                cmt = '';
        end
            
    else
        cmt = '';
    end
    
    if thinning_mode == 1; mik_mode = '_mik'; else mik_mode = ''; end
    if thinning_mode == 2; thin_mode = '_trad'; else thin_mode = ''; end
    switch tapers_mode
        case 1; long_mode = '_lt';
        case 2; long_mode = '_ltjk';
        case 3; long_mode = '_t11';
        case 4; long_mode = '_t35';
        case 5; long_mode = '_t59';
        case 7; long_mode = '_t1019';
        case 8; long_mode = '_t23';
        case 9; long_mode = '_t47';
        otherwise; long_mode = '';
    end
    %if tapers_mode == 1; long_mode = '_lt'; elseif tapers_mode == 2; long_mode = '_ltjk'; else long_mode = ''; end
    if baseline_subtract == 1; bsmode = '_bs'; else bsmode = ''; end
    if permutation_test == 1; permu_tst = '_permu'; elseif permutation_test == 2; permu_tst = '_bootstrp'; else permu_tst = ''; end
    switch coh_debias_mode
        case {1,5}
            bias_crct = '_dbsPS';    % Debias pairs of ctgs, single cell
        case {2,6}
            bias_crct = '_dbsPA';    % Debias pairs of ctgs, all cells
        case {3,7}
            bias_crct = '_dbsAS';    % Debias all ctgs, single cell
        case {4,8}
            bias_crct = '_dbsAA';    % Debias all ctgs, all ctgs cells
        otherwise
            bias_crct = '';
    end
    if coh_debias_mode > 4
        bias_crct = strcat(bias_crct,'d');      % Signifying that we're dropping trials now.
    end
    
    if do_partial_coherence == 1; ptch = '_prtl'; else ptch = ''; end
    
    switch Nwind_mode
        case 0
            nwm = '';       % Window size 300
        case 1
            nwm = '_w600';  % Window size 600
        case 2
            nwm = '_w800';  % Window size 800
        case 3
            nwm = '_w1000';  % Window size 1000
    end
    
    fname_suffix = [adj_mode cmt mik_mode thin_mode long_mode bsmode permu_tst bias_crct ptch nwm];
    
    %do_adjacent, ctgsetli_mode,  thinning_mode, tapers_mode, baseline_subtract, permutation_test
    
    
    s.do_adjacent = do_adjacent;
    s.ctgsetli_mode = ctgsetli_mode;
    s.thinning_mode = thinning_mode;
    s.tapers_mode = tapers_mode;
    s.baseline_subtract = baseline_subtract;
    s.permutation_test = permutation_test;
    s.coh_debias_mode = coh_debias_mode;
    s.do_partial_coherence = do_partial_coherence;
    
    
end


function cmt = ctgsetli_mode_2_string(ctgsetli_mode,ctgsetli_mode2)

    if nargin < 2; ctgsetli_mode2 = 0; end
    
    switch ctgsetli_mode2
        case {0,1,4}
            switch ctgsetli_mode
                case 1
                    cmt = '_sw';
                case 2
                    cmt = '_RT';     % Reaction time for match trials
                case 3
                    cmt = '_RT2';    % As RT above, but also includes match congruent and match incongruent trial groupings
                case 4
                    cmt = '_bndry';
                case 5
                    cmt = '_bndABrel';
                case 6
                    cmt = '_bndAvsB';
                case 7
                    cmt = '_bndctg60';
                case 8
                    cmt = '_nonmtch_incgt';  % Non-match incongruent trials
                case 9
                    cmt = '_bndctgfull';
                otherwise
                    cmt = '';
            end
        case {2,3}              % Second set of Ctgsetli_mode possibilities
            switch ctgsetli_mode
                case 4
                    cmt = '_bndBal';
                case 5
                    cmt = '_bndABrelBal';
                case 6
                    cmt = '_bndABcentralized';      % The mode recommended via reviewer
                case 7
                    cmt = '_bnd60ABsplit';          % 60% but split into SchA and B
                case 8
                    cmt = '_bnd50ABseparatecenter'; % As ctgsetli mode 5 but do separate analysis for center morphs
            end
    end
end