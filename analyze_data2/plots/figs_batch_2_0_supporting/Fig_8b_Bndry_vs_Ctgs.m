

function mygr = Fig_8b_Bndry_vs_Ctgs(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi,Fig8b_mode)


    plot_on = 1;
    do_sort=1;      % Sorting of data for scatterplot in Fig8b_mode = 4

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
    [sp] = plot_permute_scatter_fast(s,out.bad_any,freqband_stats,opts_PPS_temp);
    clear s

    %
    % Build group
    myN_criteria = 5;
    grt = Grp;
    grt.criteria = [2*ones(1,myN_criteria)]; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
    grt.xlims_desired = [0 120]; grt.xdata = out.group(1).xdata;


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
    
    % Sch A vs Sch B preferred
    i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryA,size(mycrit,1),1)]; mygr(i).legend = 'DA Ba';
    i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryA,size(mycrit,1),1)]; mygr(i).legend = 'DA NBa';
    i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryB,size(mycrit,1),1)]; mygr(i).legend = 'DB Bb';
    i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryB,size(mycrit,1),1)]; mygr(i).legend = 'DB NBb';
    i=i+1; mycrit = mycrit_nondecider; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryA,size(mycrit,1),1)]; mygr(i).legend = 'ND Ba';
    i=i+1; mycrit = mycrit_nondecider; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryA,size(mycrit,1),1)]; mygr(i).legend = 'ND NBa';
    i=i+1; mycrit = mycrit_nondecider; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryB,size(mycrit,1),1)]; mygr(i).legend = 'ND Bb';
    i=i+1; mycrit = mycrit_nondecider; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryB,size(mycrit,1),1)]; mygr(i).legend = 'ND NBb';
    
    % Sch A pref vs non-pref
    i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryA,size(mycrit,1),1)]; mygr(i).legend = 'DA Ba';
    i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryA,size(mycrit,1),1)]; mygr(i).legend = 'DA NBa';
    i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryA_irr,size(mycrit,1),1)]; mygr(i).legend = 'DA BaIrr';
    i=i+1; mycrit = mycrit_deciderA; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryA_irr,size(mycrit,1),1)]; mygr(i).legend = 'DA NBaIrr';
    
    % Sch B pref vs non-pref
    i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryB,size(mycrit,1),1)]; mygr(i).legend = 'DB Bb';
    i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryB,size(mycrit,1),1)]; mygr(i).legend = 'DB NBb';
    i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_bndryB_irr,size(mycrit,1),1)]; mygr(i).legend = 'DB BbIrr';
    i=i+1; mycrit = mycrit_deciderB; mygr(i) = grt; mygr(i).criteria=[mycrit]; mygr(i).ctgs = [repmat(myctgs_nonbndryB_irr,size(mycrit,1),1)]; mygr(i).legend = 'DB NBbIrr';
    
    
    for i = 1:length(mygr); [mygr(i).cells, mygr(i).numcells]= get_grouped_cells(mygr(i),[sp]); end

    [mygr] = remove_bads_cells(out.bad_any,mygr);                 % Remove bad cells from cell lists

    for i = 1:length(mygr) [mygr(i).data] = get_grouped_data(mygr(i).ctgs,mygr(i).cells,out.pls); end

    mygr = grouppairs_merge(mygr,perm2pls_dophi);

%     i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup(1:end/2)],[],opts_PM3Dcs,[]);
%     i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup(end/2+1:end)],[],opts_PM3Dcs,[]);

%     i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup(:)],[],opts_PM3Dcs,[]);
    
    if plot_on
        if Fig8b_mode == 1
            
            i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygr(1:4)],[],opts_PM3Dcs,[]);
        elseif Fig8b_mode == 2
            
            i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygr(5:6)],[],opts_PM3Dcs,[]);
            i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygr(7:8)],[],opts_PM3Dcs,[]);
        end
    end

    % 
    % groupls = out.group;
    % plsls = out.pls_stats;
    % bad_ls = out.bad_any;

    %
    
    % For description of modes, see seciton of figs_batch_2_0 entitled "Fig 8b - Ctgsetli mode 0 vs mode 5"
    if Fig8b_mode == 3
        
        % Get boundary vs non-boundary differences
        ba = out.bad_any;
        pls = out.pls;
        pls = -1*diff(pls,[],3);
        pls = pls(:,:,1:2:end);
        if perm2pls_dophi; pls = get_circ_diff(pls); end

        pls = pls * -1;
%         figure; plott_matrix3D(pls(:,~ba,:))

        ind = find_closest(out.abscissa,mean(freqband_stats));
        pls_dat = squeeze(pls(ind,:,:));



        mean(pls_dat(~ba,:))
        mygr2 = mygr(1:2);
        mygr2(1).legend= 'SchA'; 
        mygr2(2).legend= 'SchB'; 

        figure; plot(pls_dat(~ba,1),pls_dat(~ba,2),'kx','LineWidth',2);
        hold on; plot_scattergroups(pls_dat(~ba,1),pls_dat(~ba,2),mygr2,~ba); xlabel('FFC bA-NbA'); ylabel('FFC bB-NbB');
        hold on; plot_crosshair(0.5);
        windsize = 0.6;
        xlim([-windsize windsize]); ylim([-windsize windsize]);
    end
    
    if Fig8b_mode == 4
        
        % % Get boundary vs non-boundary differences
        ba = out.bad_any;
        pls = out.pls;
        pls = -1*diff(pls,[],3);
        pls = pls(:,:,1:2:end);
        if perm2pls_dophi; pls = get_circ_diff(pls); end

        pls = pls * -1;
%         figure; plott_matrix3D(pls(:,~ba,:))

        ind = find_closest(out.abscissa,mean(freqband_stats));
        pls_dat = squeeze(pls(ind,:,:));

        mean(pls_dat(~ba,:))
        
        % % Get Ctg differences (sorting)
        s = wrkspc_buffer.(['per_' mode2modename(perm_mode2)]).(stagename(curr_stage_sfc));
        Cave1 = s.Cave1; Cave2 = s.Cave2;
        
        [IA] = ind_pls(Cave1(:,:,1),Cave2(:,:,1),out.abscissa,freqband_stats);  % Sch A
        [IB] = ind_pls(Cave1(:,:,2),Cave2(:,:,2),out.abscissa,freqband_stats);  % Sch B
        [IAirr] = ind_pls(Cave1(:,:,3),Cave2(:,:,3),out.abscissa,freqband_stats);  % Sch A irr
        [IBirr] = ind_pls(Cave1(:,:,4),Cave2(:,:,4),out.abscissa,freqband_stats);  % Sch B irr
        
        pls_sort =  zeros(size(Cave1,1),size(Cave1,2),10);        
        pls_sort(:,:,1:2:end) = Cave1;
        pls_sort(:,:,2:2:end) = Cave2;
        if do_sort
            pls_sort(:,:,[1,2]) = sort_index_pls(Cave1(:,:,1),Cave2(:,:,1),IA); % Sch A Relevant
            pls_sort(:,:,[3,4]) = sort_index_pls(Cave1(:,:,2),Cave2(:,:,2),IB); % Sch B Relevant
            pls_sort(:,:,[5,6]) = sort_index_pls(Cave1(:,:,3),Cave2(:,:,3),IAirr); % Sch A IrrRelevant
            pls_sort(:,:,[7,8]) = sort_index_pls(Cave1(:,:,4),Cave2(:,:,4),IBirr); % Sch B IrrRelevant
        end
        
        % Take difference (preferred - non preferred)
        pls_sort = -1*diff(pls_sort,[],3);
        pls_sort = pls_sort(:,:,1:2:end);
        
        % Take freqband stats
        ind = find_closest(out.abscissa,mean(freqband_stats));
        pls_sortd = squeeze(pls_sort(ind,:,:));
        
        % Setup groups
        mygr2 = mygr(1:2);
        mygr2(1).legend= 'SchA'; 
        mygr2(2).legend= 'SchB'; 
        
        % Plot the differences
        
        
        % Sch A ctgs vs Sch A Bndry
        figure; % subplot(221); 
        plott_fit(pls_sortd(~ba,1),pls_dat(~ba,1),'k.');
        hold on; plot_scattergroups(pls_sortd(~ba,1),pls_dat(~ba,1),mygr2,~ba); xlabel('FFC Ctg1-Ctg2'); ylabel('FFC bA-NbA');
        
        % Sch A ctgs vs Sch B Bndry
        figure; %subplot(222); 
        plott_fit(pls_sortd(~ba,1),pls_dat(~ba,2),'k.');
        hold on; plot_scattergroups(pls_sortd(~ba,1),pls_dat(~ba,2),mygr2,~ba); xlabel('FFC Ctg1-Ctg2'); ylabel('FFC bB-NbB');
        
        % Sch A ctgs vs Sch B Bndry
        figure; %subplot(223); 
        plott_fit(pls_sortd(~ba,2),pls_dat(~ba,1),'k.');
        hold on; plot_scattergroups(pls_sortd(~ba,2),pls_dat(~ba,1),mygr2,~ba); xlabel('FFC Ctg3-Ctg4'); ylabel('FFC bA-NbA');
        
        % Sch A ctgs vs Sch B Bndry
        figure; %subplot(224); 
        plott_fit(pls_sortd(~ba,2),pls_dat(~ba,2),'k.');
        hold on; plot_scattergroups(pls_sortd(~ba,2),pls_dat(~ba,2),mygr2,~ba); xlabel('FFC Ctg3-Ctg4'); ylabel('FFC bB-NbB');

        
    end
    
end


function plot_crosshair(sz)
    if nargin < 1
        sz = 1;
    end
    plot ([-sz sz; 0 0]',[0 0; -sz sz]','k');
end