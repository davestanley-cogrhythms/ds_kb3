function [wrkspc_buffer, mygroup, pls_sort, bad_any ] = Fig_9b_Royfig_redo(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode_sort,perm_mode_sp,curr_stage_sp,freqband_stats_perm,opts_perm,do_unitsCI_vs_electBNDRY,plot_on_func)
   


    cell_group_mode = 3;    % 2-Grouped deciders; 3-all cells (as 3 for Fig9c - discrepency fixed)
    do_AB_symmetry = 1;
    sort_on = 1; % Setting to 0 completely disables all sorting
        do_sort_irr = 1; % 1-process irrelevants separately; 0-sort everything based on relevant category(maybe biased)
    
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


    % Group stuff
    mycrit = [2*ones(1,size(sp,2))];
    grt = group0(1);
    grt.criteria = []; grt.criteria_alt = []; grt.criteria_sfc = [];
    grt.xdata = group0(1).xdata;
    grt.data_name = [];

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

        case 3

            gall = group0;

            mycrit = [2*ones(1,size(sp,2))];
            for i = 1:length(gall)
                gall(i).criteria=[mycrit];
            end
        case 2
            if ~do_unitsCI_vs_electBNDRY
                % % % % % % Original, normal mode
                % Sch A preferred
                for i = 1:7; gall(i) = grt; gall(i).criteria=[mycritA]; gall(i).ctgs = i; end

                % Sch B preferred
                for i = 8:14; gall(i) = grt; gall(i).criteria=[mycritB]; gall(i).ctgs = i; end

                % Sch A Irr (i.e. Scheme A categorizations when Scheme B is cued - i.e. scheme A irrelevant)
                for i = 15:21; gall(i) = grt; gall(i).criteria=[mycritAirr]; gall(i).ctgs = i; end

                % Sch B Irr (i.e. Scheme B categorizations when Scheme A is cued - i.e. scheme B irrelevant)
                for i = 22:28; gall(i) = grt; gall(i).criteria=[mycritBirr]; gall(i).ctgs = i; end

                % Non-preferreds
                gall_np = gall(1:14);
                for i = 1:7; gall_np(i).criteria = mycritB; % Sch B NP
                end

                for i = 8:14; gall_np(i).criteria = mycritA; % Sch A NP
                end
                gall_np = gall_np([8:14,1:7]);  % Swap so Sch A non-preferred comes first

                % Non-preferrds become 29-35 (SchA) and 36-42 (SchB)
                gall = [gall gall_np];

            else
                % % % % % % Units Vs LFP boundary trials mode
                % Sch A preferred
                for i = 1:7; gall(i) = grt; gall(i).criteria=[mycritA]; gall(i).ctgs = i; end

                % Sch B preferred
                for i = 8:14; gall(i) = grt; gall(i).criteria=[mycritB]; gall(i).ctgs = i; end

                % Sch A Irr (i.e. Scheme A categorizations when Scheme B is cued - i.e. scheme A irrelevant)
                for i = 15:21; gall(i) = grt; gall(i).criteria=[mycritAirr]; gall(i).ctgs = i; end

                % Sch B Irr (i.e. Scheme B categorizations when Scheme A is cued - i.e. scheme B irrelevant)
                for i = 22:28; gall(i) = grt; gall(i).criteria=[mycritBirr]; gall(i).ctgs = i; end

                % Non-preferreds
                gall_np = gall(1:14);
                for i = 1:7; gall_np(i).criteria = non_decidersA;
                end

                for i = 8:14; gall_np(i).criteria = non_decidersB;
                end

                gall = [gall gall_np];

            end

    end
    
%     % If doing perm mode based on ctgsetli mode 7 (60% and 100%) make sure
%     % we select based on 100% morphs
%     [~, mode_subgroups] = decode_sfc_mode(perm_mode_sp);
%     [fname_suffix, do_adjacent, ctgsetli_mode, thinning_mode, tapers_mode, baseline_subtract, permutation_test, ue_pairs, coh_debias_mode, do_partial_coherence, ctgsetli_mode2] = build_sfcmode(sfc_mode, mode_subgroups);
%     if ctgsetli_mode == 7 && any(ctgsetli_mode2 == [0,1])
%         swap_map = [1:4; 5:8 ];             % FOR CTGSETLI MODE 7
% 
%         gall = group_swap.crit_move (gall,swap_map);      % THIS BUGGERS UP WHEN USING CTG 11
%         %gall = group_swap.ctgs (gall,swap_map);
% 
%     end

    for i = 1:length(gall); gall(i).legend = mylegend{mod(i-1,length(mylegend))+1}; end

    mytitles={'A Preferred Rel','B Preferred Rel','A Pref Irrel','B Pref Irrel','A Non-pref Rel','B Non-pref Rel'};
    
    if do_unitsCI_vs_electBNDRY
        mytitles = {'A Bnd A Rel','B Bnd B Rel','A Bnd A Irrel','B Bnd B Irrel','A Non-Bnd A Rel','B Non-Bnd B Rel'};
    end
    mygroup = gall;


    %
    if do_AB_symmetry
        del_ps = 7;
        clear gsymm
        for i = 1:size(pls_sort,3)/4
            gsymm(i) = group_merge(gall([i,i+del_ps]));          % Merge 1..7 and 8..14  - Preferred Relevant
            gsymm(i+7) = group_merge(gall([i+2*del_ps,i+3*del_ps])); % Merge 15...21 and 22...28 - Preferred irrelevant
            if cell_group_mode == 2 % Non-preferred relevant case is only when we have cells with a preference!
                gsymm(i+14) = group_merge(gall([i+4*del_ps,i+5*del_ps])); % Merge 29...35 and 37...42 - Non-preferred relevant (if applicable)
            end
        end
        
        %Build legend & titles
        if cell_group_mode == 2         % Grouped cells
            for i = 1:length(gsymm); gsymm(i).legend = mylegend{mod(i-1,length(mylegend))+1}; end
            mytitles={'Pref Rel','Pref Irrel','Non-pref Rel'};
        elseif cell_group_mode == 3     % All cells
            for i = 1:length(gsymm); gsymm(i).legend = mylegend{mod(i-1,length(mylegend))+1}; end
            mytitles={'All Rel','All Irrel'};
        end
        mygroup = gsymm;

        
        if do_unitsCI_vs_electBNDRY
            mytitles={'Bnd Rel','Bnd Irrel','Non-Bnd Rel'};
        end
    end



    mygroup = get_grouped_cells_and_data(mygroup,sp,pls_sort,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D);


    % Plot

    if plot_on_func
        
        % Estimate ylims
        try
            % For later versions of matlab
            mudatastats = arrayfun(@(s) mean(s.datastats),mygroup);
        catch
            % For earlier versions
            mudatastats = mygroup.myfunc_out(@(s) mean(s.datastats));
        end
        
        mylims = [min(mudatastats)*0.9 max(mudatastats)*1.1];

        % Xlims
        myxlims = [0 8];

        % Do plot

        opts_PSC.hmask = false(7);
        opts_PSC.show_legend = true;
        if do_AB_symmetry
            i=0;
            i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(1:7)],opts_PSC); ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(1).numcells))]);  xlim(myxlims);
            i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(8:14)],opts_PSC); ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(8).numcells))]); xlim(myxlims);
            if length(mygroup) > 14     % We're doing ctg deciders (cell group mode == 2)
                i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(15:21)],opts_PSC); ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(15).numcells))]); xlim(myxlims);
            end

        else
            i=0;
            i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(1:7)],opts_PSC);ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(1).numcells))]);  xlim(myxlims);
            i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(8:14)],opts_PSC); ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(8).numcells))]);  xlim(myxlims);
            i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(15:21)],opts_PSC); ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(15).numcells))]);  xlim(myxlims);
            i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(22:28)],opts_PSC); ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(22).numcells))]);  xlim(myxlims);
            if cell_group_mode == 2
                i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(29:35)],opts_PSC); ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(29).numcells))]);  xlim(myxlims);
                i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(36:42)],opts_PSC); ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(36).numcells))]);  xlim(myxlims);
            end
        end
    end
end
