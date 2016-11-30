

function mygr = Fig_8c_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)


    plot_on = 1;
    plot_on_stats = 1;

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
    opts_PPS_temp.split_plusminus = 0;

    % Get statistics
    % Do sp for perm_mode (ctg sensitivity)
    s = wrkspc_buffer.(['per_' mode2modename(perm_mode2)]).(stagename(curr_stage_sfc));
    sp = plot_permute_scatter_fast(s,out.bad_any,freqband_stats,opts_PPS_temp);
    clear s

    %
    % Build group
    myN_criteria = 5;
    grt = Grp;
    grt.criteria = [2*ones(1,myN_criteria)]; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
    grt.xlims_desired = [0 120]; grt.xdata = out.group(1).xdata;
    grt.interpreter = 'none';


                % Criteria
                mycrit_deciderA = [1 0 2 2 2;];   % Sch A and/or B deciders
                mycrit_deciderB = [0 1 2 2 2;];   % Sch A and/or B deciders
                mycrit_nondecider = [0 0 2 2 2];

                % Relevant Ctgs
                myctgs_bndryA = [2];         % Sch A Rel 50 % trials
                myctgs_nonbndryA = [1];      % Sch A Rel 100 % trials
                myctgs_bndryB = [4];         % Sch B Rel 50 % trials
                myctgs_nonbndryB = [3];      % Sch B Rel 100 % trials

                % Irrelevant Ctgs
                myctgs_bndryA_irr = [8];         % IrrRel Sch A 50 % trials
                myctgs_nonbndryA_irr = [7];      % IrrRel Sch A 100 % trials
                myctgs_bndryB_irr = [6];         % IrrRel Sch B 50 % trials
                myctgs_nonbndryB_irr = [5];      % IrrRel Sch B 100 % trials

            %      % New, implemented August 12
            %         ctgsetli = [ [t100 & schA] [t50 & schA] [t100 & schB] [t50 & schB] ...
            %                 [t100_irr & schA] [t50_irr & schA] [t100_irr & schB] [t50_irr & schB] ...
            %             ];

        
    clear mygroup
    i=0;
    
    % Sch A Deciders; SchA Bndry Rel vs SchB Bndry Rel vs SchB Bndry Irrel
    i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryA,size(mycrit,1),1)]; mygr(i).legend = 'dA bA';
    i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryA,size(mycrit,1),1)]; mygr(i).legend = 'dA nbA';
    i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryB,size(mycrit,1),1)]; mygr(i).legend = 'dA bB';
    i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryB,size(mycrit,1),1)]; mygr(i).legend = 'dA nbB';
    i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryB_irr,size(mycrit,1),1)]; mygr(i).legend = 'dA bBirr';
    i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryB_irr,size(mycrit,1),1)]; mygr(i).legend = 'dA nbBirr';
    
    % Sch B Deciders; SchB Bndry Rel vs SchA Bndry Rel vs SchA Bndry Irrel
    i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryB,size(mycrit,1),1)]; mygr(i).legend = 'dB bB';
    i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryB,size(mycrit,1),1)]; mygr(i).legend = 'dB nbB';
    i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryA,size(mycrit,1),1)]; mygr(i).legend = 'dB bA';
    i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryA,size(mycrit,1),1)]; mygr(i).legend = 'dB nbA';
    i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryA_irr,size(mycrit,1),1)]; mygr(i).legend = 'dB bAirr';
    i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryA_irr,size(mycrit,1),1)]; mygr(i).legend = 'dB nbAirr';
    
    % Nondeciders; SchA Bndry Rel vs SchB Bndry Rel vs SchB Bndry Irrel
    i=i+1; mycrit = mycrit_nondecider; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryA,size(mycrit,1),1)]; mygr(i).legend = 'NP bA';
    i=i+1; mycrit = mycrit_nondecider; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryA,size(mycrit,1),1)]; mygr(i).legend = 'NP nbA';
    i=i+1; mycrit = mycrit_nondecider; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryB,size(mycrit,1),1)]; mygr(i).legend = 'NP bB';
    i=i+1; mycrit = mycrit_nondecider; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryB,size(mycrit,1),1)]; mygr(i).legend = 'NP nbB';
    i=i+1; mycrit = mycrit_nondecider; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryB_irr,size(mycrit,1),1)]; mygr(i).legend = 'NP bBirr';
    i=i+1; mycrit = mycrit_nondecider; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryB_irr,size(mycrit,1),1)]; mygr(i).legend = 'NP nbBirr';
    
    
    for i = 1:length(mygr); [mygr(i).cells, mygr(i).numcells]= get_grouped_cells(mygr(i),[sp]); end

    [mygr] = remove_bads_cells(out.bad_any,mygr);                 % Remove bad cells from cell lists

    for i = 1:length(mygr)
        [mygr(i).data] = get_grouped_data(mygr(i).ctgs,mygr(i).cells,out.pls);
        [mygr(i).metadata.mypairs] = get_grouped_data(1,mygr(i).cells,out.mypairs')';
        [mygr(i).datastats, mygr(i).freqband_stats] = calc_pls_stats(out.abscissa,mygr(i).data,freqband_stats,'do_mean_ctgs',1);
        
    end

    mygr = grouppairs_merge(mygr,perm2pls_dophi);
    
    mygr(1).legend = 'Cat/Dog Bndry';   % Pref A
    mygr(2).legend = 'Fat/Thin Bndry';
    mygr(4).legend = 'Fat/Thin Bndry';  % Pref B
    mygr(5).legend = 'Cat/Dog Bndry';
    mygr(7).legend = 'Cat/Dog Bndry';  % Nonpref
    mygr(8).legend = 'Thin/Fat Bndry';
    
    
    %'$\rm \bf {Pref_{A}} Ctg 1v2$';

%     i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup(1:end/2)],[],opts_PM3Dcs,[]);
%     i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup(end/2+1:end)],[],opts_PM3Dcs,[]);

%     i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup(:)],[],opts_PM3Dcs,[]);
    
    if plot_on
        i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygr([1:2])],[],opts_PM3Dcs,[]);
        i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygr([4:5])],[],opts_PM3Dcs,[]);
        % Do non pref as well
        %i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygr([1:2,7:8])],[],opts_PM3Dcs,[]);
        %i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygr([4:5,7:8])],[],opts_PM3Dcs,[]);
    end
    
    if plot_on_stats
        opts_PSC.hmask = blkdiag(true(2),true(2));
        opts_PSC.remove_dependent=0;
        opts_PSC.paperfig_mode=1;
        i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(mygr([1,2,4,5]),opts_PSC);
    end

    % 
    % groupls = out.group;
    % plsls = out.pls_stats;
    % bad_ls = out.bad_any;

    %
    
 
    
    
end