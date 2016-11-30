

function mygroup = Fig_8a_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)

    plot_on = 1;

    % Plot matrix settings
    opts_PM3Dcs.paperfig_mode = 1;
    opts_PM3Dcs.do_sgolay = 0;
    opts_PM3Dcs.stats_mode = 0;

    % Statistics settings
    opts_PPS_temp.paperfig_mode = 1;
    opts_PPS_temp.hide_boot_results = 1;
    opts_PPS_temp.hide_all_plots = 1;
    opts_PPS_temp.do_bh0=1;
    opts_PPS_temp.do_phi=perm2pls_dophi;
    opts_PPS_temp.split_plusminus=0;

    % Get statistics
    % Do sp for perm_mode (ctg sensitivity)
    s = wrkspc_buffer.(['per_' mode2modename(perm_mode2)]).(stagename(curr_stage_sfc));
    [sp] = plot_permute_scatter_fast(s,out.bad_any,freqband_stats,opts_PPS_temp);
    clear s

    %
    % Build group
    myN_criteria = 5;
    gr_template = Grp;
    gr_template.criteria = [2*ones(1,myN_criteria)]; gr_template.criteria_alt = []; gr_template.criteria_sfc = []; gr_template.ctgs = 1;
    gr_template.xlims_desired = [0 120]; gr_template.xdata = out.group(1).xdata;


    % Criteria
    mycrit_deciderA = [1 0 2 2 2; 0 1 2 2 2; 1 1 2 2 2;];   % Sch A and/or B deciders
    mycrit_nondecider = [0 0 2 2 2];
    
    
    
    % Boundary / Non-boundary Ctgs
    myctgs_bndryA = [6];         % Sch A Rel 50 % trials
    myctgs_nonbndryA = [5];      % Sch A Rel 100 % trials


    
        
%         ctgsetli = [[ctgsetli_subset & t100_mat] [ctgsetli_subset & t60_mat] ...    % 40% & 60% morphs rel
%             [ctgsetli_subset & t100_irr_mat] [ctgsetli_subset & t60_irr_mat] ...    % 40% & 60% morphs irrel
%             [ctgsetli_subset & t100_mat] [ctg_correct_or_noresponse & t50_mat] ...            % 50% morphs rel
%             [ctgsetli_subset & t100_irr_mat] [ctg_correct_or_noresponse & t50_irr_mat] ...    % 50% morphs irrrel
%             [ctgsetli_subset & t100_mat] [ctgsetli_subset & near_bnd] ...            % 50% morphs rel
%             [ctgsetli_subset & t100_irr_mat] [ctgsetli_subset & near_bnd_irr] ...    % 50% morphs irrrel
%             ];  

    clear group
    i=0;
    
    % Sch A vs Sch B preferred
    
    i=i+1; mycrit = mycrit_deciderA; mygroup(i) = gr_template; mygroup(i).criteria=[mycrit]; mygroup(i).ctgs = [repmat(myctgs_bndryA,size(mycrit,1),1)]; mygroup(i).legend = 'Pref';
    i=i+1; mycrit = mycrit_deciderA; mygroup(i) = gr_template; mygroup(i).criteria=[mycrit]; mygroup(i).ctgs = [repmat(myctgs_nonbndryA,size(mycrit,1),1)]; mygroup(i).legend = '';
    
    i=i+1; mycrit = mycrit_nondecider; mygroup(i) = gr_template; mygroup(i).criteria=[mycrit]; mygroup(i).ctgs = [repmat(myctgs_bndryA,size(mycrit,1),1)]; mygroup(i).legend = 'No-pref';
    i=i+1; mycrit = mycrit_nondecider; mygroup(i) = gr_template; mygroup(i).criteria=[mycrit]; mygroup(i).ctgs = [repmat(myctgs_nonbndryA,size(mycrit,1),1)]; mygroup(i).legend = '';
    
    
    
    
    for i = 1:length(mygroup); [mygroup(i).cells, mygroup(i).numcells]= get_grouped_cells(mygroup(i),[sp]); end

    [mygroup] = remove_bads_cells(out.bad_any,mygroup);                 % Remove bad cells from cell lists

    for i = 1:length(mygroup)
        [mygroup(i).data] = get_grouped_data(mygroup(i).ctgs,mygroup(i).cells,out.pls);
        [mygroup(i).metadata.mypairs] = get_grouped_data(ones(size(mygroup(i).cells,2),1),mygroup(i).cells,out.mypairs',out.funames1D)';
        [mygroup(i).datastats, mygroup(i).freqband_stats] = calc_pls_stats(out.abscissa,mygroup(i).data,freqband_stats,'do_mean_ctgs',1);
    end
    

    mygroup = grouppairs_merge(mygroup,perm2pls_dophi);
    
    mygroup(1).legend = 'Categorizers';
    mygroup(2).legend = 'Non-categorizers';
    
    if plot_on

    %     i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup(1:end/2)],[],opts_PM3Dcs,[]);
    %     i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup(end/2+1:end)],[],opts_PM3Dcs,[]);

    %     i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup(:)],[],opts_PM3Dcs,[]);

        i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup],[],opts_PM3Dcs,[]);
        
        opts_PSC.hmask = blkdiag(true(2));
        i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup],opts_PSC);
    end

    % 
    % groupls = out.group;
    % plsls = out.pls_stats;
    % bad_ls = out.bad_any;

    %
    clear out
end