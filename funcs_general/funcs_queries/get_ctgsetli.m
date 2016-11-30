
function [ctgsetli,N_ctgs_base,N_ctgs_extras] = get_ctgsetli(sfc_mode,md,ctgsetli_mode,ctgsetli_mode2, filenum)
    %% [ctgsetli,N_ctgs_base,N_ctgs_extras] = get_ctgsetli(sfc_mode,md,ctgsetli_mode)
    
    
    % Set up some stuff that's used universally
    sotr_match = md.tr_match(md.sample_on_trials); sotr_nonmatch = md.tr_nonmatch(md.sample_on_trials);
    sot = md.sample_on_trials;
    etr = md.each_trial_response;
    
    %ctg_all = ones(size(ctgsetli,1),1);
    ctg_all = get_good_samples(1:length(md.sample_on_trials),md);   % I am updating ctg_all to only include correct trials. Prior to analyze_data7, this also included incorrect trials.
    ctg_all = ctg_all(:);
    ctg_temp = etr == 0; ctg_correct = ctg_temp(sot);     % Should equal ctg_all
    ctg_temp = etr == 0 | etr == 1; ctg_correct_or_noresponse = ctg_temp(sot);
                    % Good or no response (no response is ok for boundary
                    % trials since they could be "non-match"; not sure why
                    % they're not marked as correct.  early is not OK!
                    % (because they could press BEFORE the test image
                    % appears).
    
    
    ctg_incorrect = md.each_trial_response == 6;
    ctg_incorrect = ctg_incorrect(md.sample_on_trials);
    
    % ******************
    % warning('correct treatment of correct vs incorrect trials... ');
    % The relevant ones should include "no response" case since nonmatch
    % doesn't rule out no response errors. Might also want to include early
    % ones. For Irrel, should only use "correct" only since correct works
    % properly here. Maybe rename ctg_all to "ctg_good" or something
    % Maybe email Jefferson.
    % ******************
    
    
    % Can think of ctgsetli_mode2 as the "tens" digit in ctgsetli_mode
    % Thus, ctgsetli_mode can range from 0 to 99. Full
    % mode is ctgsetli_mode2*10 + ctgsetli_mode
    % 0-9   - normal options
    % 10-19 - as 0-9, but combine rel and irrel
    
    
    switch ctgsetli_mode2
        case {0,1,4}
                  % Modes 0-9 and 10-19
                  % 10-19 are 0-9 plus Mrg

            if ctgsetli_mode == 0 || ctgsetli_mode == 1 || ctgsetli_mode == 7 || ctgsetli_mode == 9
                %% Ctgsetli normal - 00, 01, 07, 09
                if floor(sfc_mode) == 5
                    N_ctgs_extras = 2;       % Test Cat/Dog; Test Match/Nonmatch; Scheme A/B
                else
                    N_ctgs_extras = 2;
                end



                if isfield(md,'ctg1_nr_trials')     % Is Roy data.
                    N_ctgs_base = 4*2;
                    N_categories = N_ctgs_base + N_ctgs_extras + 1;
                    ctgsetli = false(length(get_good_samples(md.ctg1_trials,md)),N_categories);
                    ctgsetli(:,1) = get_good_samples(md.ctg1_trials,md);
                    ctgsetli(:,2) = get_good_samples(md.ctg2_trials,md);
                    ctgsetli(:,3) = get_good_samples(md.ctg3_trials,md);
                    ctgsetli(:,4) = get_good_samples(md.ctg4_trials,md);

                    ctgsetli(:,5) = get_good_samples(md.ctg1_nr_trials,md);
                    ctgsetli(:,6) = get_good_samples(md.ctg2_nr_trials,md);
                    ctgsetli(:,7) = get_good_samples(md.ctg3_nr_trials,md);
                    ctgsetli(:,8) = get_good_samples(md.ctg4_nr_trials,md);
                else
                    N_ctgs_base = 4*2;
                    N_categories = N_ctgs_base + N_ctgs_extras + 1;
                    ctgsetli = false(length(get_good_samples(md.ctg1_trials,md)),N_categories);
                    ctgsetli(:,1) = get_good_samples(md.ctg1_trials,md);
                    ctgsetli(:,2) = get_good_samples(md.ctg2_trials,md);
                    ctgsetli(:,3) = get_good_samples(md.ctg3_trials,md);
                    ctgsetli(:,4) = get_good_samples(md.ctg4_trials,md); 


                    % Dummy data for cromer
                    ctgsetli(:,5:8) = false(size(ctgsetli,1),4);
                    ctgsetli(1:5,[5,7]) = true;
                    ctgsetli(6:10,[6,8]) = true;

                end

                %if floor(sfc_mode) == 5
                if 1
                    ctg1test = (ctgsetli(:,1) & sotr_match) | (ctgsetli(:,2) & sotr_nonmatch);
                    ctg2test = (ctgsetli(:,2) & sotr_match) | (ctgsetli(:,1) & sotr_nonmatch);
                    ctg3test = (ctgsetli(:,3) & sotr_match) | (ctgsetli(:,4) & sotr_nonmatch);
                    ctg4test = (ctgsetli(:,4) & sotr_match) | (ctgsetli(:,3) & sotr_nonmatch);
                    schA_trials = ctgsetli(:,1) | ctgsetli(:,2);
                    schB_trials = ctgsetli(:,3) | ctgsetli(:,4);
                    %ctgsetli_extras = [ctg1test, ctg2test, ctg3test, ctg4test, sotr_match, sotr_nonmatch, schA_trials, schB_trials];
                    ctgsetli_extras = [schA_trials, schB_trials];
                else
                    ctgsetli_extras = [];
                end

                % Format base
                ctgsetli_base = ctgsetli(:,1:8);

                % Collapse base if combining rel vs irrel
                if ctgsetli_mode2 == 1
                    % Combine rel and irrel
                    ctgsetli_base(:,1:4) = [any(ctgsetli_base(:,[1,5]),2), any(ctgsetli_base(:,[2,6]),2), ...
                                    any(ctgsetli_base(:,[3,7]),2), any(ctgsetli_base(:,[4,8]),2)];            
                    ctgsetli_base(:,5:8) = false(size(ctgsetli,1),4);
                end

%%
                if ctgsetli_mode == 0
                    %% Ctgsetli 00 - Combine everything
                    ctgsetli = [ctgsetli_base, ctgsetli_extras, ctg_all];

                    if ctgsetli_mode2 == 1  % If merge rel/irrel, add dummy data
                        %% Ctgsetli 10
                        d1 = false(size(ctgsetli,1),1); d2=d1;
                        d1(1:5,:) = true; d2(6:10,:) = true;
                        ctgsetli(:,5:8) = [d1,d2,d1,d2];
                    end
                    
                    if ctgsetli_mode2 == 4  % Swap things so that rel and irrel are paired
                        %% Ctgsetli 40
                        ctgsetli(:,1:8) = ctgsetli(:,[1,5,2,6,3,7,4,8]);
                    end

                    N_ctgs_base = size(ctgsetli_base,2);
                    N_ctgs_extras = size(ctgsetli_extras,2);

                elseif ctgsetli_mode == 1
                    %% Ctgsetli 01 - Modify for switch trials
                    N_ctgs_base = 6;
                    N_ctgs_extras = 0;
                    tsw = get_switch_trials(md);
                    tsw=tsw(:);
                    %ctgsetli_subset = ctgsetli(:,9:11);
                    %tsw_mat = repmat(tsw(:),1,size(ctgsetli_subset,2));
                    ctgsetli = [[ctgsetli_extras(:,1) & ~tsw] [ctgsetli_extras(:,1) & tsw] ...
                        [ctgsetli_extras(:,2) & ~tsw] [ctgsetli_extras(:,2) & tsw] ...
                        [ctg_all & ~tsw] [ctg_all & tsw] ...
                        ];   % 1-3 SchA/B/all & non-switch; 4-6 SchA/B/all & switch

                elseif ctgsetli_mode == 7
                    %% Ctgsetli 07 - Boundry trials for Ctg1 vs Ctg2
                        % Note, this seems to be highly redundant with ctgsetli_mode 9.
                        % Can perhaps delete. - actually, it's still useful because
                        % we can do stats on this!
                    N_ctgs_base = 16;
                    N_ctgs_extras = 0;
                    [t60, t60_irr, t50, t50_irr, t80, t80_irr, t100, t100_irr] = get_boundary_trials(md);

                    sot = md.sample_on_trials;
                    etr = md.each_trial_response;

                    % Matify
                    t60_mat = repmat(t60(:),1,size(ctgsetli_base,2));
                    t100_mat = repmat(t100(:),1,size(ctgsetli_base,2));


                    ctgsetli = [ [t60_mat & ctgsetli_base] ...
                                [t100_mat & ctgsetli_base] ...
                                ];

                    if ctgsetli_mode2 == 1  % If merge rel/irrel, add dummy data
                        %% Ctgsetli 17
                        d1 = false(size(ctgsetli,1),1); d2=d1;
                        d1(1:5,:) = true; d2(6:10,:) = true;
                        ctgsetli(:,5:8) = [d1,d2,d1,d2];
                        ctgsetli(:,13:16) = [d1,d2,d1,d2];
                    end

                elseif ctgsetli_mode == 9
                    %% %% Ctgsetli 09 - Ctg1 vs Ctg2 for various boundary conditions - highly redundant with mode 7
                    N_ctgs_base = 28;
                    N_ctgs_extras = 0;
                    [t60, t60_irr, t50, t50_irr, t80, t80_irr, t100, t100_irr] = get_boundary_trials(md);

                    sch = md.sch;
                    schA = strcmp(sch,'A');
                    schB = strcmp(sch,'B');

                    sot = md.sample_on_trials;

                    % Shorten to sample on
                    schA = schA(sot)';      
                    schB = schB(sot)';


                    % Correct SchA and schB, force to be either "correct" or "no response"
                    %    This rules out other bad trials, such as broken
                    %    fixation
                    ctgsetli_extras2 = [schA(:) & ctg_correct_or_noresponse schB(:) & ctg_correct_or_noresponse];   % Correction, made Aug 17th

                    % Matify
                    t60_mat = repmat(t60(:),1,size(ctgsetli_base,2));
                    t80_mat = repmat(t80(:),1,size(ctgsetli_base,2));
                    t100_mat = repmat(t100(:),1,size(ctgsetli_base,2));
                    t50_mat = repmat(t50(:),1,size(ctgsetli_extras2,2));
                    t50_mat_irr = repmat(t50_irr(:),1,size(ctgsetli_extras2,2));


                    ctgsetli = [ [t60_mat & ctgsetli_base] ...
                                 [t80_mat & ctgsetli_base] ...
                                 [t100_mat & ctgsetli_base] ...
                                 [t50_mat & ctgsetli_extras2] ...
                                 [t50_mat_irr & ctgsetli_extras2] ...
                                ];

                    N_ctgs_base = size(ctgsetli,2);

                    if ctgsetli_mode2 == 1  % If merge rel/irrel, add dummy data
                        %% Ctgsetli 19 
                        d1 = false(size(ctgsetli,1),1); d2=d1;
                        d1(1:5,:) = true; d2(6:10,:) = true;
                        ctgsetli(:,5:8) = [d1,d2,d1,d2];
                        ctgsetli(:,13:16) = [d1,d2,d1,d2];
                        ctgsetli(:,21:24) = [d1,d2,d1,d2];
                    end

                end


            elseif ctgsetli_mode == 2
                %% %% Ctgsetli 02 - Short and long reaction time trials
                N_ctgs_base = 2;
                N_ctgs_extras = 0;
                N_categories = N_ctgs_base + N_ctgs_extras + 1;
                plot_on = 0;

                % Extract reaction times for quantile processing
                    % Note - this is ONLY used for getting quantile thresholds
                rt0 = md.match_rt; 
                rt0 = rt0 + mean(rt0(~isnan(rt0)))*1e-4*randn(size(rt0)); % Add a small amount of noise to reaction times.

                rt = rt0(sotr_match &  ctg_all);
                temp=sum(isnan(rt));    % Should be zero NaNs
                if temp >0
                    warning('Isnans presnet in reaction time data');
                end

                % Get quantiles
                p0 = [0.25 0.5 0.75];
                % p0 = [0.333, 0.666];        
                q = calc_quantiles(p0,rt);

                % Get trials associated with short and long RTs
                ctg_rt_short = ctg_all & sotr_match & rt0 < q(1);   % Good trials + match trials + short RT
                ctg_rt_long = ctg_all & sotr_match & rt0 > q(end);  % Good trials + match trials + long RT

                rt_short = sort(md.match_rt(ctg_rt_short));
                rt_long = sort(md.match_rt(ctg_rt_long));
                number_of_nans = sum(isnan(rt_short)) + sum(isnan(rt_long));
                if number_of_nans > 0
                    warning('NaNs detected in RT trials.');
                end

                if plot_on
                    rt_all = sort(md.match_rt);
                    figure; plot(rt_all);
                    hold on; plot(rt_short,'r');
                    hold on; plot(rt_long,'g');
                    legend('All RTs','Short RTs','Long RTs');
                    title(['Number of NaNs= ' num2str(number_of_nans) ]);
                end

                ctgsetli = [ctg_rt_short, ctg_rt_long, ctg_all];

            elseif ctgsetli_mode == 3
                %% %% Ctgsetli 03 - Short and long reaction time trials with congruent vs incongruent
                N_ctgs_base = 4;
                N_ctgs_extras = 0;
                plot_on = 0;

                condspath = getpath('path_conditions');
                load(fullfile(condspath,'conditions_organized.mat'))

        %         match_congruent = match_congruent{filenum};
        %         match_incongruent = match_incongruent{filenum};

                mc = get_good_samples(match_congruent{filenum},md);
                mic = get_good_samples(match_incongruent{filenum},md);

                % Extract reaction times for quantile processing - congruent
                clear rt rt0
                mc_curr = mc;
                rt0 = md.match_rt;
                rt0 = rt0 + mean(rt0(~isnan(rt0)))*1e-4*randn(size(rt0)); % Add a small amount of noise to reaction times.
                                                            % Without this,
                                                            % if there are repeated values, these will all be lumped together.
                                                            % This will create uneven sample sizes. With this noise present, it
                                                            % will cause a subset of these repeated values to be randomly
                                                            % selected, enabling equal sample size.
                rt = rt0(mc_curr' &  ctg_all);
                temp=sum(isnan(rt));    % Should be zero NaNs
                if temp >0
                    warning('Isnans presnet in reaction time data');
                end

                % Get quantiles
                p0 = [0.25 0.5 0.75];
                % p0 = [0.333, 0.666];        
                q = calc_quantiles(p0,rt);


                % Get trials associated with short and long RTs
                mc_short = ctg_all & mc_curr' & rt0 < q(1);   % Good trials + match trials + short RT
                mc_long = ctg_all & mc_curr' & rt0 > q(end);  % Good trials + match trials + long RT

                rt_short = sort(md.match_rt(mc_short));
                rt_long = sort(md.match_rt(mc_long));
                number_of_nans = sum(isnan(rt_short)) + sum(isnan(rt_long));
                if number_of_nans > 0
                    warning('NaNs detected in RT trials.');
                end


                % Extract reaction times for quantile processing - incongruent
                clear rt rt0
                mc_curr = mic;
                rt0 = md.match_rt;
                rt0 = rt0 + mean(rt0(~isnan(rt0)))*1e-4*randn(size(rt0)); % Add a small amount of noise to reaction times.

                rt = rt0(mc_curr' &  ctg_all);
                temp=sum(isnan(rt));    % Should be zero NaNs
                if temp >0
                    warning('Isnans presnet in reaction time data');
                end

                % Get quantiles
                p0 = [0.25 0.5 0.75];
                % p0 = [0.333, 0.666];        
                q = calc_quantiles(p0,rt);


                % Get trials associated with short and long RTs
                mic_short = ctg_all & mc_curr' & rt0 < q(1);   % Good trials + match trials + short RT
                mic_long = ctg_all & mc_curr' & rt0 > q(end);  % Good trials + match trials + long RT

                rt_short = sort(md.match_rt(mic_short));
                rt_long = sort(md.match_rt(mic_long));
                number_of_nans = sum(isnan(rt_short)) + sum(isnan(rt_long));
                if number_of_nans > 0
                    warning('NaNs detected in RT trials.');
                end


                if plot_on
                    rt_all = sort(md.match_rt);
                    figure; plot(rt_all);
                    hold on; plot(rt_short,'r');
                    hold on; plot(rt_long,'g');
                    legend('All RTs','Short RTs','Long RTs');
                    title(['Number of NaNs= ' num2str(number_of_nans) ]);
                end

                ctgsetli = [mc_short, mc_long, mic_short, mic_long];


            elseif ctgsetli_mode == 4
                %% Ctgsetli 04 - For boundary trials
                N_ctgs_base = 12;
                N_ctgs_extras = 0;
                [t60, t60_irr, t50, t50_irr, t80, t80_irr, t100, t100_irr] = get_boundary_trials(md);

                t60 = t60(:) & ctg_correct;
                t60_irr = t60_irr(:) & ctg_correct;
                t50 = t50(:) & ctg_correct_or_noresponse;
                t50_irr = t50_irr(:) & ctg_correct_or_noresponse;
                t100 = t100(:) & ctg_correct;
                t100_irr = t100_irr(:) & ctg_correct;


                near_bnd = t60 | t50;
                near_bnd_irr = t60_irr | t50_irr;

                t50_allowed = t50 & t100_irr;
                t50_irr_allowed = t50_irr & t100;
                t60_complement = t100_irr | t60_irr;
                t60_irr_complement = t100 | t60;



                ctgsetli = [[t100] [t60] ...    % 40% & 60% morphs rel
                    [t100_irr] [t60_irr] ...    % 40% & 60% morphs irrel
                    [t100] [ctg_correct_or_noresponse & t50] ...            % 50% morphs rel
                    [t100_irr] [ctg_correct_or_noresponse & t50_irr] ...    % 50% morphs irrrel
                    [t100] [near_bnd] ...            % 50% morphs rel
                    [t100_irr] [near_bnd_irr] ...    % 50% morphs irrrel
                    ];  


                % We consider no response trials correct for the 50% cases also. We
                % consider this because (using the test_boundary_codes function
                % below) we ascertained that non-match trials are NEVER classified
                % as no response. Therefore, since ALL boundary trials could
                % potentially be non-matches, it is not right to EVER discount no
                % response as being correct!

                if ctgsetli_mode2 == 1
                    %% Ctgsetli 14
                    % Doesn't really make sense to combine relevant and
                    % irrelevant boundary trials, but here we go!
                    warning('Probably not proper to do this. 100% and 50% can potentially both contain the same trial due to combination.');

                    % Add dummy data
                    d1 = false(size(ctgsetli,1),1); d2=d1;
                    d1(1:5,:) = true; d2(6:10,:) = true;

                    ctgsetli = [any(ctgsetli(:,[1,3]),2), any(ctgsetli(:,[2,4]),2), ...
                                d1, d2, ...
                                any(ctgsetli(:,[5,7]),2), any(ctgsetli(:,[6,8]),2), ...
                                d1, d2, ...
                                any(ctgsetli(:,[9,11]),2), any(ctgsetli(:,[10,12]),2), ...
                                d1, d2];
                end


            elseif ctgsetli_mode == 5
                %% Ctgsetli 05 - Boundry trials for SchA and SchB
                N_ctgs_base = 12;
                N_ctgs_extras = 0;
                [t60, t60_irr, t50, t50_irr, t80, t80_irr, t100, t100_irr] = get_boundary_trials(md);


                sch = md.sch;
                morphA = md.morphA;
                morphB = md.morphB;

                schA = strcmp(sch,'A');
                schB = strcmp(sch,'B');

                sot = md.sample_on_trials;
                etr = md.each_trial_response;

                % Shorten to sample on
                etr = etr(sot);
                schA = schA(sot)';      
                schB = schB(sot)';

                % Only good.
                t60 = t60(:) & ctg_correct;
                t60_irr = t60_irr(:) & ctg_correct;
                t50 = t50(:) & ctg_correct_or_noresponse;
                t50_irr = t50_irr(:) & ctg_correct_or_noresponse;
                t100 = t100(:) & ctg_correct;
                t100_irr = t100_irr(:) & ctg_correct;

                tcentral = t60 | t50;
                tcentral_irr = t60_irr | t50_irr;

        %         ctgsetli = [ [t100 & schA] [t60 & schA] [t100 & schB] [t60 & schB] ...
        %                 [t100 & schA] [t50 & schA] [t100 & schB] [t50 & schB] ...
        %                 [t100 & schA] [tcentral & schA] [t100 & schB] [tcentral & schB] ...
        %             ];

                % New, implemented August 12
                ctgsetli = [ [t100 & schA] [t50 & schA] [t100 & schB] [t50 & schB] ...
                        [t100_irr & schA] [t50_irr & schA] [t100_irr & schB] [t50_irr & schB] ...
                    ];

                if ctgsetli_mode2 == 1      
                    %% Ctgsetli 15
                            % Doesn't really make sense to combine relevant and
                            % irrelevant boundary trials, but here we go!
                    % Add dummy data
                    d1 = false(size(ctgsetli,1),1); d2=d1;
                    d1(1:5,:) = true; d2(6:10,:) = true;
                    
                    ctgsetli = [any(ctgsetli(:,[1,7]),2), any(ctgsetli(:,[2,8]),2), any(ctgsetli(:,[3,5]),2), any(ctgsetli(:,[4,6]),2) ...
                            d1,d2,d1,d2];
                end



            elseif ctgsetli_mode == 6
                %% Ctgsetli 06 - Boundry trials for SchA vs SchB
                N_ctgs_base = 6;
                N_ctgs_extras = 0;
                [t60, t60_irr, t50, t50_irr, t80, t80_irr, t100, t100_irr] = get_boundary_trials(md);

                sch = md.sch;
                morphA = md.morphA;
                morphB = md.morphB;

                schA = strcmp(sch,'A');
                schB = strcmp(sch,'B');

                sot = md.sample_on_trials;
                etr = md.each_trial_response;

                % Shorten to sample on
                etr = etr(sot);
                schA = schA(sot)';      
                schB = schB(sot)';

                % Only good.
                t60 = t60(:) & ctg_correct;
                t60_irr = t60_irr(:) & ctg_correct;
                t50 = t50(:) & ctg_correct_or_noresponse;
                t50_irr = t50_irr(:) & ctg_correct_or_noresponse;
                t100 = t100(:) & ctg_correct;
                t100_irr = t100_irr(:) & ctg_correct;

                tcentral = t60 | t50;
                tcentral_irr = t60_irr | t50_irr;

                ctgsetli = [ [t60 & schA] [t60 & schB] ...
                           [t50 & schA] [t50 & schB] ...
                           [tcentral & schA] [tcentral & schB] ...
                    ];

                if ctgsetli_mode2 == 1
                    % We're only considering rel, so not applicable here. In
                    % fact, I don't think we really use this mode anyways!
                end


            elseif ctgsetli_mode == 8
                %% Ctgsetli 08 - Nonmatch + congruent vs incongruent
                N_ctgs_base = 2;
                N_ctgs_extras = 0;

                condspath = getpath('path_conditions');
                load(fullfile(condspath,'conditions_organized.mat'))

                nmc = get_good_samples(nonmatch_congruent{filenum},md);
                nmic = get_good_samples(nonmatch_incongruent{filenum},md);



                ctgsetli = [ [nmc(:) & ctg_all] [nmic(:) & ctg_all]  ...
                    ];

                if ctgsetli_mode2 == 1
                    % We're only considering rel, so not applicable here. In
                    % fact, I don't think we really use this mode anyways!
                end
            end
        case {2,3}
            switch ctgsetli_mode
                case 4
                    %% Ctgsetli 24 - For boundary trials
                    N_ctgs_base = 12;
                    N_ctgs_extras = 0;

                    sch = md.sch;
                    morphA = md.morphA;
                    morphB = md.morphB;

                    schA = strcmp(sch,'A');
                    schB = strcmp(sch,'B');

                    sot = md.sample_on_trials;

                    % Shorten to sample on
                    schA = schA(sot)';      
                    schB = schB(sot)';
                    morphA = morphA(sot);
                    morphB = morphB(sot);



                    [t60, t60_irr, t50, t50_irr, t80, t80_irr, t100, t100_irr] = get_boundary_trials(md);

                    t60 = t60(:) & ctg_correct;
                    t60_irr = t60_irr(:) & ctg_correct;
                    t50 = t50(:) & ctg_correct_or_noresponse;
                    t50_irr = t50_irr(:) & ctg_correct_or_noresponse;
                    t100 = t100(:) & ctg_correct;
                    t100_irr = t100_irr(:) & ctg_correct;


                    near_bnd = t60 | t50;
                    near_bnd_irr = t60_irr | t50_irr;

                    t50_allowed = t50 & t100_irr;
                    t50_irr_allowed = t50_irr & t100;
                    t60_complement = t100_irr | t60_irr;
                    t60_irr_complement = t100 | t60;

    %                 % old code
    %                 ctgsetli = [[t100] [t60] ...    % 40% & 60% morphs rel
    %                     [t100_irr] [t60_irr] ...    % 40% & 60% morphs irrel
    %                     [t100] [ctg_correct_or_noresponse & t50] ...            % 50% morphs rel
    %                     [t100_irr] [ctg_correct_or_noresponse & t50_irr] ...    % 50% morphs irrrel
    %                     [t100] [near_bnd] ...            % 50% morphs rel
    %                     [t100_irr] [near_bnd_irr] ...    % 50% morphs irrrel
    %                     ];  

                    ctgsetli = [[t100 & t60_complement] [t60 & t60_complement] ...    % 40% & 60% morphs rel
                        [t100_irr & t60_irr_complement] [t60_irr & t60_irr_complement] ...    % 40% & 60% morphs irrel
                        [t100 & t100_irr] [ctg_correct_or_noresponse & t50_allowed] ...            % 50% morphs rel
                        [t100_irr & t100] [ctg_correct_or_noresponse & t50_irr_allowed] ...    % 50% morphs irrrel
                        ];  

                    ctgsetli = [ctgsetli, ...               % ctgsetli(:,9..12) - 50% or 60%
                        [ctgsetli(:,1) | ctgsetli(:,5)], [ctgsetli(:,2) | ctgsetli(:,6)] ...         % 100% vs 60% or 50%
                        [ctgsetli(:,3) | ctgsetli(:,7)], [ctgsetli(:,4) | ctgsetli(:,8)] ...         % 100% irr vs 60% or 50% irr
                        ];


                    % We consider no response trials correct for the 50% cases also. We
                    % consider this because (using the test_boundary_codes function
                    % below) we ascertained that non-match trials are NEVER classified
                    % as no response. Therefore, since ALL boundary trials could
                    % potentially be non-matches, it is not right to EVER discount no
                    % response as being correct!

                    if ctgsetli_mode2 == 3
                        %% Ctgsetli 34
                        
                        % Add dummy data
                        d1 = false(size(ctgsetli,1),1); d2=d1;
                        d1(1:5,:) = true; d2(6:10,:) = true;

                        ctgsetli = [any(ctgsetli(:,[1,3]),2), any(ctgsetli(:,[2,4]),2), ...
                                    d1, d2, ...
                                    any(ctgsetli(:,[5,7]),2), any(ctgsetli(:,[6,8]),2), ...
                                    d1, d2, ...
                                    any(ctgsetli(:,[9,11]),2), any(ctgsetli(:,[10,12]),2), ...
                                    d1, d2];
                    end


                case 5
                    %% Ctgsetli 25 - Boundry trials for SchA and SchB
                    N_ctgs_base = 8;
                    N_ctgs_extras = 0;
                    [t60, t60_irr, t50, t50_irr, t80, t80_irr, t100, t100_irr] = get_boundary_trials(md);


                    sch = md.sch;
                    morphA = md.morphA;
                    morphB = md.morphB;

                    schA = strcmp(sch,'A');
                    schB = strcmp(sch,'B');

                    sot = md.sample_on_trials;
                    etr = md.each_trial_response;

                    % Shorten to sample on
                    etr = etr(sot);
                    schA = schA(sot)';      
                    schB = schB(sot)';
                    morphA = morphA(sot);
                    morphB = morphB(sot);

                    % Only good.
                    t60 = t60(:) & ctg_correct;
                    t60_irr = t60_irr(:) & ctg_correct;
                    t50 = t50(:) & ctg_correct_or_noresponse;
                    t50_irr = t50_irr(:) & ctg_correct_or_noresponse;
                    t100 = t100(:) & ctg_correct;
                    t100_irr = t100_irr(:) & ctg_correct;

    %                 % Old method, not caring about symmetry
    %                 ctgsetli = [ [t100 & schA] [t50 & schA] ...
    %                              [t100 & schB] [t50 & schB] ...
    %                              [t100_irr & schA] [t50_irr & schA] ...
    %                              [t100_irr & schB] [t50_irr & schB] ...
    %                     ];


                    ctgsetli = [ [t100 & schA & (t100_irr | t50_irr)] [t50 & schA] ...
                                 [t100 & schB & (t100_irr | t50_irr)] [t50 & schB] ...
                                 [t100_irr & schA & (t100 | t50)] [t50_irr & schA] ...
                                 [t100_irr & schB & (t100 | t50)] [t50_irr & schB] ...
                        ];

                    if ctgsetli_mode2 == 3      % Combine relevant and irrelevant
                        %% Ctgsetli 35
                        % Add dummy data
                        d1 = false(size(ctgsetli,1),1); d2=d1;
                        d1(1:5,:) = true; d2(6:10,:) = true;

                        ctgsetli = [any(ctgsetli(:,[1,7]),2), any(ctgsetli(:,[2,8]),2), any(ctgsetli(:,[3,5]),2), any(ctgsetli(:,[4,6]),2) ...
                            d1,d2,d1,d2];
                    end
                    
                case 6
                    %% Ctgsetli 26 - Boundry trials for SchA and SchB
                    N_ctgs_base = 12;
                    N_ctgs_extras = 0;
                    [t60, t60_irr, t50, t50_irr, t80, t80_irr, t100, t100_irr] = get_boundary_trials(md);


                    sch = md.sch;
                    morphA = md.morphA;
                    morphB = md.morphB;

                    schA = strcmp(sch,'A');
                    schB = strcmp(sch,'B');

                    sot = md.sample_on_trials;
                    etr = md.each_trial_response;

                    % Shorten to sample on
                    etr = etr(sot);
                    schA = schA(sot)';      
                    schB = schB(sot)';
                    morphA = morphA(sot);
                    morphB = morphB(sot);

                    % Only good.
                    t60 = t60(:) & ctg_correct;
                    t60_irr = t60_irr(:) & ctg_correct;
                    t50 = t50(:) & ctg_correct_or_noresponse;
                    t50_irr = t50_irr(:) & ctg_correct_or_noresponse;
                    t100 = t100(:) & ctg_correct;
                    t100_irr = t100_irr(:) & ctg_correct;

                    ctgsetli = [ [t50 & t50_irr & schA] [t50 & ~t50_irr & schA] ...
                                 [t50 & t50_irr & schB] [t50 & ~t50_irr & schB] ...
                                 [t50 & t50_irr & schA] [~t50 & t50_irr & schA] ...
                                 [t50 & t50_irr & schB] [~t50 & t50_irr & schB] ...
                                 [t50 & ~t50_irr & schA] [t100 & t100_irr & schA] ...
                                 [t50 & ~t50_irr & schB] [t100 & t100_irr & schB] ...
                               ];

                    if ctgsetli_mode2 == 3      % Combine relevant and irrelevant
                        %% Ctgsetli 36
                        
                        d1 = false(size(ctgsetli,1),1); d2=d1;
                        d1(1:80,:) = true; d2(81:160,:) = true;

                        
                        ctgsetli = [any(ctgsetli(:,[1,7]),2), any(ctgsetli(:,[2,8]),2), any(ctgsetli(:,[3,5]),2), any(ctgsetli(:,[4,6]),2) ...
                            d1,d2,d1,d2, ...
                            any(ctgsetli(:,[9,11]),2), any(ctgsetli(:,[10,12]),2), ...
                            d1,d2 ...
                            ];
                    end
                case 7
                    %% Ctgsetli 27 - For 60% morphs - version of ctgsetli mode 5 for 60%
                    N_ctgs_base = 8;
                    N_ctgs_extras = 0;

                    sch = md.sch;
                    morphA = md.morphA;
                    morphB = md.morphB;

                    schA = strcmp(sch,'A');
                    schB = strcmp(sch,'B');

                    sot = md.sample_on_trials;

                    % Shorten to sample on
                    schA = schA(sot)';      
                    schB = schB(sot)';
                    morphA = morphA(sot);
                    morphB = morphB(sot);

                    [t60, t60_irr, t50, t50_irr, t80, t80_irr, t100, t100_irr] = get_boundary_trials(md);

                    t60 = t60(:) & ctg_correct;
                    t60_irr = t60_irr(:) & ctg_correct;
                    t100 = t100(:) & ctg_correct;
                    t100_irr = t100_irr(:) & ctg_correct;

                    t60_complement = t100_irr | t60_irr;
                    t60_irr_complement = t100 | t60;
                    
                    ctgsetli = [[t60 & t60_complement & schA] [t100 & t60_complement & schA] ...
                                [t60 & t60_complement & schB] [t100 & t60_complement & schB] ...
                                [t60_irr & t60_irr_complement & schA] [t100_irr & t60_irr_complement & schA] ...
                                [t60_irr & t60_irr_complement & schB] [t100_irr & t60_irr_complement & schB] ...
                               ];
                           
                    if ctgsetli_mode2 == 3      % Combine relevant and irrelevant
                        %% Ctgsetli 37
                        
                        d1 = false(size(ctgsetli,1),1); d2=d1;
                        d1(1:5,:) = true; d2(6:10,:) = true;

                        
                        ctgsetli = [any(ctgsetli(:,[1,7]),2), any(ctgsetli(:,[2,8]),2), any(ctgsetli(:,[3,5]),2), any(ctgsetli(:,[4,6]),2) ...
                            d1,d2,d1,d2, ...
                            ];
                    end
                    
                case 8
                    %% Ctgsetli 28 - For 50% morphs - as ctgsetli mode 5, but do seprate analysis for centre
                    N_ctgs_base = 12;
                    N_ctgs_extras = 0;

                    sch = md.sch;
                    morphA = md.morphA;
                    morphB = md.morphB;

                    schA = strcmp(sch,'A');
                    schB = strcmp(sch,'B');

                    sot = md.sample_on_trials;

                    % Shorten to sample on
                    schA = schA(sot)';      
                    schB = schB(sot)';
                    morphA = morphA(sot);
                    morphB = morphB(sot);

                    [t60, t60_irr, t50, t50_irr, t80, t80_irr, t100, t100_irr] = get_boundary_trials(md);

                    t50 = t50(:) & ctg_correct_or_noresponse;
                    t50_irr = t50_irr(:) & ctg_correct_or_noresponse;
                    t100 = t100(:) & ctg_correct;
                    t100_irr = t100_irr(:) & ctg_correct;

                    t50_allowed = t50 & t100_irr;
                    t50_irr_allowed = t50_irr & t100;

                    ctgsetli = [[t50 & t100_irr & schA] [t100 & t100_irr & schA] ...     % SchA Bnd vs Outer
                                [t50 & t100_irr & schB] [t100 & t100_irr & schB] ...     % SchB Bnd vs Outer
                                [t50_irr & t100 & schA] [t100 & t100_irr & schA] ...     % SchA Bnd vs Outer Irr
                                [t50_irr & t100 & schB] [t100 & t100_irr & schB] ...     % SchB Bnd vs Outer Irr
                                [t50 & t50_irr & schA] [t100 & t100_irr & schA] ...     % SchA Center vs Outer
                                [t50 & t50_irr & schB] [t100 & t100_irr & schB] ...     % SchB Center vs Outer
                                ];
                                

                    if ctgsetli_mode2 == 3      % Combine relevant and irrelevant
                        %% Ctgsetli 38
                        
                        d1 = false(size(ctgsetli,1),1); d2=d1;
                        d1(1:5,:) = true; d2(6:10,:) = true;

                        
                        ctgsetli = [any(ctgsetli(:,[1,7]),2), any(ctgsetli(:,[2,8]),2), any(ctgsetli(:,[3,5]),2), any(ctgsetli(:,[4,6]),2) ...
                            d1,d2,d1,d2, ...
                            any(ctgsetli(:,[9,11]),2), any(ctgsetli(:,[10,12]),2), d1, d2, ...
                            ];
                    end

                    

            end
    end
end


function q = calc_quantiles(p0,rt)

    plot_on = 0;
    
    q = zeros(1,length(p0));
    for i = 1:length(p0)
        p_curr=p0(i);
        rt = sort(rt);
        q(i) = quantile(rt,p_curr);
        if plot_on;
            %% Test plot
            figure; [f,x]=ecdf(rt); plot(x,f); hold on; plot([q(i),q(i)],[0,1],'r');
            legend('empirical CDF',['Quantile @ ' num2str(p_curr)]); ylabel('probability'); xlabel('Reaction Time');
        end
    end
end





function test_boundary_codes(md,t50,t50_irr)
    %% Test plot - test boundary codes
    etr = md.each_trial_response;
    sot = md.sample_on_trials;
    m = md.tr_match;
    nm = md.tr_nonmatch;

    etr = etr(sot);
    m = m(sot);
    nm = nm(sot);

    figl;
    subplot(411); hist(etr(t50),0:10); title('Boundary');
    subplot(412); hist(etr(t50_irr),0:10); title('Boundary Irr');
    subplot(413); hist(etr(m),0:10); title('Match trials');
    subplot(414); hist(etr(nm),0:10); title('Nonmatch trials');


end