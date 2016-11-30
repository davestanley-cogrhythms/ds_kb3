

function [ctgsetli_A, ctgsetli_good_A] = calc_actual_good_ctgsetli(sfc_mode,bad_trials)
    % Estimates the ctgsetli values for all cells & electrode pairs,
    % including the estimated good trials
    
    if ~exist('sfc_mode','var'); sfc_mode = 22.4013110; end
    if ~exist('bad_trials','var'); 
        stage=3;    % (For testing.._
        bad_trials = prep_badtrials(stage);
    end

    %% Setup MD, variables
    exclude_clipping_trials = 1;
    
    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix, do_adjacent, ctgsetli_mode, thinning_mode, tapers_mode, baseline_subtract, permutation_test, ue_pairs, coh_debias_mode, do_partial_coherence, ctgsetli_mode2] = build_sfcmode(sfc_mode, mode_subgroups);
    sfc_mode_group = floor(sfc_mode); 
    
    path_md = fullfile(getpath('path_metadata'),'sansunits');
    file_list = get_filelist(path_md);
    md_all = cellfun(@(f) load (fullfile(path_md,f)),file_list);     % Load metadata
    
    
    
    %% Loop through all files, and get ctgsetli, ctgsetli_good
    Nfiles = length(md_all);
%     tic
    for i = 1:Nfiles
        %i
        
        md = md_all(i);
        filenum = i;
        
        Nunits = length(md.unit_names);    % Number of units tracked
        Nelects = length(md.lfp_names);      % Number of electrodes
        


        % Set up mypairs for iteration
        switch ue_pairs
            case {0,2};                                             % Current unit versus current/adjacent field (SFC)
                Ncoherences = Nunits;                           
                mypairs = [1:(Nunits); 1:(Nunits)]';
            case 3                                                  % Current unit all fields (SFC)
                % Generate all possible unit-electrode pairs
                [units1,lfp2] = meshgrid(1:Nunits,1:Nelects);
                mypairs = [units1(:) lfp2(:)];
                %mypairs = mypairs(units1 ~= 0,:);
                clear units1 lfp2
                Ncoherences = size(mypairs,1);
            case 4                                                  % Current field versus all field (FFC)
                mypairs = upper_triangle_pairs(Nelects);
                Ncoherences = size(mypairs,1);
            case 5                                                  % Current unit versus all units (SSC)
                mypairs = upper_triangle_pairs(Nunits);
                Ncoherences = size(mypairs,1);
            case 6
                mypairs = [1:(Nelects); 1:(Nelects)]';
                Ncoherences = size(mypairs,1);
        end
        
        [ctgsetli0] = get_ctgsetli(sfc_mode,md,ctgsetli_mode,ctgsetli_mode2,filenum);
        
        for j = 1:Ncoherences
            curr_unit = j;
            %j;

            % Identify clipping trials
            if exclude_clipping_trials
                switch ue_pairs
                    case 0
                        curr_lfp = md.unit_LFP_pairs(2,curr_unit);
                        curr_badtrials = bad_trials{filenum}(:,curr_lfp);
                    case 2
                        [curr_lfp]= get_adjacent_electrode(md,curr_unit);
                        curr_badtrials = bad_trials{filenum}(:,curr_lfp);
                    case 3
                        curr_lfp = mypairs(j,2);
                        curr_badtrials = bad_trials{filenum}(:,curr_lfp);
                    case 4
                        curr_lfp = mypairs(j,1);
                        curr_lfp2 = mypairs(j,2);
                        curr_badtrials = bad_trials{filenum}(:,curr_lfp) | bad_trials{filenum}(:,curr_lfp2);
                    case 5
                        curr_badtrials = false(Nsamples,1); % All trials are good!
                    case 6
                        curr_lfp = j;
                        curr_badtrials = bad_trials{filenum}(:,curr_lfp);
                end
            else
                curr_badtrials = false(Nsamples,1); % All trials are good!
            end

            ctgsetli=ctgsetli0;
            if i==29 && (ue_pairs ~= 5)
                %fprintf(['Shortened file #29: ' md.recording_name '. Should match L112106.mat.\n']);
                ctgsetli = ctgsetli(1:length(curr_badtrials),:);       % Shorten ctgsetli if we're missing lfp data. Should only be for file #29. Only do this if we are consider LFP
            end
            
            
            if ue_pairs ~= 5
                ctgsetli_good = ctgsetli & repmat(~curr_badtrials(:),1,size(ctgsetli,2));
            else
                ctgsetli_good = ctgsetli;   % If doing unit-unit comparison, there is no clipping and all trials are good
            end
                
            
            ctgsetli_A{i}{j} = sum(ctgsetli);
            ctgsetli_good_A{i}{j} = sum(ctgsetli_good);
            
        end
                
    end
%     toc
    


end