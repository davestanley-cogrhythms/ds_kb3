function [wrkspc_buffer, gall, gall2 ] = Fig_9c_Cat_Index(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode_sort,perm_mode_sp,curr_stage_sp,freqband_stats_perm,opts_perm,do_unitsCI_vs_electBNDRY,plot_on_func,do_AB_symmetry)

    plot_regression = 0;

    cell_group_mode = 3;    % 1-enumerate ctgs (obsolete); 2-Grouped deciders; 3-All cells (same as 1 for Fig9b)
    %do_AB_symmetry = 1;
    sort_on = 0; % Setting to 0 completely disables all sorting
        do_sort_irr = 1; % 1-process irrelevants separately; 0-sort everything based on relevant category(maybe biased)
        
    do_group_collapse_pls2days = 0;
    
    if do_unitsCI_vs_electBNDRY
        cell_group_mode = 2;        % Force to group mode 2
    end
    
    [wrkspc_buffer, out_pls, out_sort, out_perm, sp] = Fig_9_general_load_and_arrange_CI(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode_sort,perm_mode_sp,curr_stage_sp,freqband_stats_perm,opts_perm,do_unitsCI_vs_electBNDRY);
    pls = out_pls.pls;
    pls_stats = out_pls.pls_stats;
    abscissa = out_pls.abscissa;
    bad_any = out_pls.bad_any;
    funames1D = out_pls.funames1D;
    mypairs = out_pls.mypairs;
    group0 = out_pls.group;
    
    pls_sort = get_pls_sort(out_pls, out_sort, sort_on, do_sort_irr, freqband_stats);

    mylegend = {'100','80','60','50','40','20','0'};
    
    del_p=8;
    del_ps = 7;
    % Convert to cat index
    pls_CI_A = morphs2CI(pls_sort(:,:,[1:7]+del_ps*0),bad_any);
    pls_CI_B = morphs2CI(pls_sort(:,:,[1:7]+del_ps*1),bad_any);
    pls_CI_A_irr = morphs2CI(pls_sort(:,:,[1:7]+del_ps*2),bad_any);
    pls_CI_B_irr = morphs2CI(pls_sort(:,:,[1:7]+del_ps*3),bad_any);
    
    pls_CI = cat(3,pls_CI_A,pls_CI_B,pls_CI_A_irr,pls_CI_B_irr);
    clear pls_CI_A pls_CI_B pls_CI_A_irr pls_CI_B_irr
    
    % Group stuff
    mycrit = [2*ones(1,size(sp,2))];
    grt = group0(1);
    grt.criteria = []; grt.criteria_alt = []; grt.criteria_sfc = [];
    grt.xdata = group0(1).xdata;
    grt.data_name = ['Category Index'];

    
    % Set up criteria definitions
    mycrit = [2*ones(1,size(sp,2))];
    mycritA = [2*ones(1,size(sp,2))];
    mycritB = [2*ones(1,size(sp,2))];
    mycritAirr = [2*ones(1,size(sp,2))];
    mycritBirr = [2*ones(1,size(sp,2))];
    if do_unitsCI_vs_electBNDRY
        mycritA(1:2) = [1 2];
        mycritB(1:2) = [2 1];
    else
        mycritA(1:2) = [1 0];
        mycritB(1:2) = [0 1];
    end
    if do_sort_irr
        mycritAirr(3:4) = [1 0];
        mycritBirr(3:4) = [0 1];
    else
        mycritAirr = mycritA;
        mycritBirr = mycritB;
    end
    non_deciders = [2*ones(1,size(sp,2))]; non_decidersA = non_deciders; non_decidersB = non_deciders;
    non_decidersA(1:2) = [0,2];
    non_decidersB(1:2) = [2,0];
    
    clear gall
    switch cell_group_mode
        case 1
            N_arrangements = 2;
            
            for i = 1:size(pls_CI,3)
                gall(i).criteria=[mycrit]; gall(i).ctgs = i;
            end
        case 2
            N_arrangements= 3;
            if ~do_unitsCI_vs_electBNDRY
                i=0;
                % Sch A preferred
                i=i+1; gall(i) = grt; gall(i).criteria=[mycritA]; gall(i).ctgs = 1; gall(i).legend = 'SchA Pref';

                % Sch B preferred
                i=i+1; gall(i) = grt; gall(i).criteria=[mycritB]; gall(i).ctgs = 2; gall(i).legend = 'SchB Pref';

                % Non-preferreds
                % Sch A non-preferred rel
                i=i+1; gall(i) = grt; gall(i).criteria=[mycritA]; gall(i).ctgs = 2; gall(i).legend = 'SchA Non-Pref';

                % Sch B non-preferred rel
                i=i+1; gall(i) = grt; gall(i).criteria=[mycritB]; gall(i).ctgs = 1; gall(i).legend = 'SchB Non-Pref';

                % Preferred irrelevant
                % Sch A preferred irr
                i=i+1; gall(i) = grt; gall(i).criteria=[mycritAirr]; gall(i).ctgs = 3; gall(i).legend = 'SchA Pref Irr';

                % Sch B preferred irr
                i=i+1; gall(i) = grt; gall(i).criteria=[mycritBirr]; gall(i).ctgs = 4; gall(i).legend = 'SchB Pref Irr';
            else
                i=0;
                % Sch A preferred
                i=i+1; gall(i) = grt; gall(i).criteria=[mycritA]; gall(i).ctgs = 1; gall(i).legend = 'A Bnd A Rel';

                % Sch B preferred
                i=i+1; gall(i) = grt; gall(i).criteria=[mycritB]; gall(i).ctgs = 2; gall(i).legend = 'B Bnd B Rel';

                % Non-preferreds
                % Sch A non-preferred rel
                i=i+1; gall(i) = grt; gall(i).criteria=[non_decidersA]; gall(i).ctgs = 1; gall(i).legend = 'A Non-Bnd A Rel';

                % Sch B non-preferred rel
                i=i+1; gall(i) = grt; gall(i).criteria=[non_decidersB]; gall(i).ctgs = 2; gall(i).legend = 'B Non-Bnd B Rel';

                % Preferred irrelevant
                % Sch A preferred irr
                i=i+1; gall(i) = grt; gall(i).criteria=[mycritAirr]; gall(i).ctgs = 3; gall(i).legend = 'A Bnd A Irr';

                % Sch B preferred irr
                i=i+1; gall(i) = grt; gall(i).criteria=[mycritBirr]; gall(i).ctgs = 4; gall(i).legend = 'B Bnd B Irr';

            end

        case 3
            
            mycrit = [2*ones(1,size(sp,2))];
            
            N_arrangements= 2;
            
            i=0;
            % Sch A preferred
            i=i+1; gall(i) = grt; gall(i).criteria=[mycrit]; gall(i).ctgs = i; gall(i).legend = 'C/D ';
            
            % Sch B preferred
            i=i+1; gall(i) = grt; gall(i).criteria=[mycrit]; gall(i).ctgs = i; gall(i).legend = 'G/D';
            
            % Sch A preferred irr
            i=i+1; gall(i) = grt; gall(i).criteria=[mycrit]; gall(i).ctgs = i; gall(i).legend = 'C/D Irr';
            
            % Sch B preferred irr
            i=i+1; gall(i) = grt; gall(i).criteria=[mycrit]; gall(i).ctgs = i; gall(i).legend = 'G/D Irr';
    end

    % If doing perm mode based on ctgsetli mode 7 (60% and 100%) make sure
    % we select based on 100% morphs
    [~, mode_subgroups] = decode_sfc_mode(perm_mode_sp);
    [fname_suffix, do_adjacent, ctgsetli_mode, thinning_mode, tapers_mode, baseline_subtract, permutation_test, ue_pairs, coh_debias_mode, do_partial_coherence, ctgsetli_mode2] = build_sfcmode(sfc_mode, mode_subgroups);
    if ctgsetli_mode == 7 && any(ctgsetli_mode2 == [0,1])
        swap_map = [1:4; 5:8 ];             % FOR CTGSETLI MODE 7

        gall = group_swap.crit_move (gall,swap_map);      % THIS BUGGERS UP WHEN USING CTG 11
        %gall = group_swap.ctgs (gall,swap_map);

    end

    %
    if do_AB_symmetry
        clear gsymm
        i=0;
        i=i+1; gsymm(i) = group_merge(gall(1),gall(2));
        i=i+1; gsymm(i) = group_merge(gall(3),gall(4));
        if cell_group_mode == 2
            i=i+1; gsymm(i) = group_merge(gall(5),gall(6));
        end
        

        %Build legend & titles
        if cell_group_mode == 2         % Grouped cells
            mytitles={'Pref Rel','Non-pref Rel','Pref Irrel'};
        elseif cell_group_mode == 3
            mytitles={'All Rel','All Irrel'};
        end
        for i = 1:length(gsymm); gsymm(i).legend = mytitles{i}; end
        
    
        if do_unitsCI_vs_electBNDRY
            %Build legend & titles
            mytitles={'Bnd Rel.','Non-Bnd Rel.','Bnd Irrel.'};
            for i = 1:length(gsymm); gsymm(i).legend = mytitles{i}; end
        end
        
    
        gall = gsymm;

    end



    gall = get_grouped_cells_and_data(gall,sp,pls_CI,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D);
    

    if do_group_collapse_pls2days
        gall = group_collapse_pls2days(gall);
    end


    % Plot
    opts_PSC.remove_dependent = 0;
    opts_PSC.show_legend = true;
    opts_PSC.hmask = [];    % Default (Just compare adjacent entries)
    

    % Produce gall2 for plotting - Rearranged version of gall 
    if ~do_unitsCI_vs_electBNDRY
        if ~do_AB_symmetry
            myordering = 2;          
            if cell_group_mode == 2
                if myordering == 1
                    gall2 = gall([1,3,5,2,4,6]); % 1-[SchA Rel, SchA Npref, SchA Irrel; SchB likewise];
                    opts_PSC.hmask = blkdiag(true(N_arrangements),true(N_arrangements));
                elseif myordering == 2
                    gall2 = gall([1,5,2,6,3,4]); % 2-[SchA Rel, SchA Irrel; SchB Rel, SchB Irrel; SchA Npref, SchB Npref]
                end
            elseif cell_group_mode == 3
                gall2 = gall([1,3,2,4]); % [SchA Rel, SchA Irrel, SchB Rel, SchB Irrel]
            end
        else    % Do symmetry
            opts_PSC.hmask = true(N_arrangements);
            gall2 = gall;
        end
    else                    % We're in do_unitsCI_vs_electBNDRY mode
        if ~do_AB_symmetry
            gall2 = gall([1,3,2,4]);
        else
            opts_PSC.hmask = true(2);
            %gall2 = gall;
            gall2 = gall(1:2);
        end
    end
        
    if plot_on_func
        
        % Do plot
        i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[gall2]);
        i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(gall2,opts_PSC);

    end
    
    if plot_regression
        
        %% Regression plots
        
        mymode=1;     % 1-convert electrode pairs to units (work in units)
                      % 2-convert units to electrode pairs (work in elects)
        
        % Pull out Y statistic
        Y = (out_perm.Cave2 - out_perm.Cave1) ;
        % Y = wrkspc_buffer.per_22_4514111.stage2.zsc_Cave; % Use zscore instead of abs difference
        
        % Y is boundary response (FFC)
        index = find_closest(out_perm.abscissa,mean(freqband_stats_perm));
        Y = squeeze(Y(index,:,:,:));
        
        % X is category index (units)
        X = pls_CI;
        index = find_closest(out_pls.abscissa,mean(freqband_stats));
        X = squeeze(X(index,:,:,:));
        
        switch mymode
            case 1
                % Map FFC to units
                ulpairs = build_unit_LFP_pairs_all(wrkspc_buffer.currmd.md);
                unit2lfp = ulpairs(2,:);
                mypairsX_units = out_pls.mypairs(:,2);
                mypairsX_elects = unit2lfp(mypairsX_units);
                mypairsY_FFC = out_perm.mypairs;
                Y2 = zeros(size(X,1),size(Y,2));
                for i = 1:length(mypairsX_elects);
                    ind = find(mypairsX_elects(i) == mypairsY_FFC(:,1) | mypairsX_elects(i) == mypairsY_FFC(:,2));
                    Y2(i,:) = mean(Y(ind,:),1);
                    %Y2(i,:) = max(Y(ind,:),[],1);
                end
                Y = Y2;
                bads_regress = any(isnan([X,Y]),2);
            case 2
                % Map FFC to units
                ulpairs = build_unit_LFP_pairs_all(wrkspc_buffer.currmd.md);
                unit2lfp = ulpairs(2,:);
                mypairsX_units = out_pls.mypairs(:,2);
                mypairsY_FFC = out_perm.mypairs;
                
                % lfp2unit = zeros(1,max(mypairsY_FFC(:)));
                % lfp2unit(:) = NaN;
                % lfp2unit(unit2lfp) = 1:length(unit2lfp);
                X2 = zeros(size(Y,1),size(X,2)) * NaN;
                for i = 1:length(mypairsY_FFC)
                    ind = find(mypairsY_FFC(i,1) == unit2lfp | mypairsY_FFC(i,2) == unit2lfp); % Find units associated with each electrode pair
                    if ~isempty(ind)
                        X2(i,:) = mean(X(ind,:),1);
                    else
                    end
                end
                X = X2;
                bads_regress = any(isnan([X,Y]),2);

        end
        
        % Plot relationship
        X=X(~bads_regress,:);
        Y=Y(~bads_regress,:);
%         X=zscore(X);
%         Y=zscore(Y);
        ind=1; figure; plott_fit(X(:,ind),Y(:,ind),'.')
        ind=2; figure; plott_fit(X(:,ind),Y(:,ind),'.')
    end

end


