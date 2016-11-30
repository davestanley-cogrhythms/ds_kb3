function mygroup = Fig_6b_Ctgs_Rel_Irrel(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)
    %% Setup
    % Settings
    plot_on = 1;
    do_AB_symmetry = 1;
    get_statistics = 1;
    
    %
    % Plot matrix settings
    opts_PM3Dcs.paperfig_mode = 1;
    opts_PM3Dcs.do_sgolay = 0;
    opts_PM3Dcs.stats_mode = 0;

    % Plot bar settings
    opts_PSC.paperfig_mode = 1;
    opts_PSC.hmask = false(7);
    opts_PSC.remove_dependent = 0;
    
    % Statistics settings
    opts_perm = Opts_Perm;
    opts_perm.do_bh0=1;
    opts_perm.do_phi=perm2pls_dophi;
    
    
    

    %% Sorting into pref / non-pref
    % Sort into preferred and non-preferred
    Cave1 = wrkspc_buffer.(['per_' mode2modename(perm_mode2)]).(stagename(curr_stage_sfc)).Cave1;
    Cave2 = wrkspc_buffer.(['per_' mode2modename(perm_mode2)]).(stagename(curr_stage_sfc)).Cave2;


    [IA] = ind_pls(Cave1(:,:,1),Cave2(:,:,1),out.abscissa,freqband_stats);  % Sch A
    [IB] = ind_pls(Cave1(:,:,2),Cave2(:,:,2),out.abscissa,freqband_stats);  % Sch B
    
    if size(Cave1,3) > 3
        [IAirr] = ind_pls(Cave1(:,:,3),Cave2(:,:,3),out.abscissa,freqband_stats);  % Sch A
        [IBirr] = ind_pls(Cave1(:,:,4),Cave2(:,:,4),out.abscissa,freqband_stats);  % Sch B
        warning('FYI: Sorting Irr deciders. This may harm stats, but be more correct');
    else
        IAirr = IA;
        IBirr = IB;
    end
    
    

    sz = size(out.pls);
    pls_sort = zeros([sz(1),sz(2),8]);
    
    
    pls_sort(:,:,[1,2]) = sort_index_pls(out.pls(:,:,1),out.pls(:,:,2),IA); % Sch A Relevant
    pls_sort(:,:,[3,4]) = sort_index_pls(out.pls(:,:,3),out.pls(:,:,4),IB); % Sch A Relevant
    pls_sort(:,:,[5,6]) = sort_index_pls(out.pls(:,:,5),out.pls(:,:,6),IAirr); % Sch A Relevant
    pls_sort(:,:,[7,8]) = sort_index_pls(out.pls(:,:,7),out.pls(:,:,8),IBirr); % Sch A Relevant
    
    mylegend = {'100','80','60','50','40','20','0'};
    mytitles = {'Pref. Rel.','NPref. Rel.','Pref. Irrel'};

    %% Get statistics
    if get_statistics
        %%
        % Do sp for perm_mode (ctg sensitivity)
        s = wrkspc_buffer.(['per_' mode2modename(perm_mode2)]).(stagename(curr_stage_sfc));
        [sp1] = plot_permute_scatter_fast(s,out.bad_any,freqband_stats,opts_perm);
        clear s
        
        sp = [sp1];
    else
        Ncrit = 5;
        sp = false(size(pls_sort,2),Ncrit);
    end

    % Setup groups
    % Build group template
    mycrit = [2*ones(1,size(sp,2))];
    gr_template = Grp;
    gr_template.criteria = mycrit; gr_template.criteria_alt = []; gr_template.criteria_sfc = []; gr_template.ctgs = 1;

    clear gall

    if size(sp,2) > 3
        mycritA = [1,0,2,2,2];
        mycritB = [0,1,2,2,2];
        mycritA_irr = [2,2,1,0,2];
        mycritB_irr = [2,2,1,0,2];
        warning('FYI: Grouping Irr deciders. This may harm stats, but be more correct');
        
    else
        mycritA = [1,0,2];
        mycritB = [0,1,2];
        mycritA_irr = mycritA;
        mycritB_irr = mycritB;
    
    end

    i=0;
    i=i+1; gall(i) = gr_template; gall(i).criteria=[mycritA]; gall(i).ctgs = 1; gall(i).legend = 'Morphs >50%';
    i=i+1; gall(i) = gr_template; gall(i).criteria=[mycritA]; gall(i).ctgs = 2; gall(i).legend = 'Morphs <50%';
    i=i+1; gall(i) = gr_template; gall(i).criteria=[mycritB]; gall(i).ctgs = 3; gall(i).legend = '';
    i=i+1; gall(i) = gr_template; gall(i).criteria=[mycritB]; gall(i).ctgs = 4; gall(i).legend = '';
    
    i=i+1; gall(i) = gr_template; gall(i).criteria=[mycritA]; gall(i).ctgs = 3; gall(i).legend = 'Morphs >50%';
    i=i+1; gall(i) = gr_template; gall(i).criteria=[mycritA]; gall(i).ctgs = 4; gall(i).legend = 'Morphs <50%';
    i=i+1; gall(i) = gr_template; gall(i).criteria=[mycritB]; gall(i).ctgs = 1; gall(i).legend = '';
    i=i+1; gall(i) = gr_template; gall(i).criteria=[mycritB]; gall(i).ctgs = 2; gall(i).legend = '';
    
    i=i+1; gall(i) = gr_template; gall(i).criteria=[mycritA_irr]; gall(i).ctgs = 5; gall(i).legend = 'Morphs >50%';
    i=i+1; gall(i) = gr_template; gall(i).criteria=[mycritA_irr]; gall(i).ctgs = 6; gall(i).legend = 'Morphs <50%';
    i=i+1; gall(i) = gr_template; gall(i).criteria=[mycritB_irr]; gall(i).ctgs = 7; gall(i).legend = '';
    i=i+1; gall(i) = gr_template; gall(i).criteria=[mycritB_irr]; gall(i).ctgs = 8; gall(i).legend = '';
    
    
    mygroup = gall;
    
    %
    if do_AB_symmetry
        clear gsymm
        
        i=0;
        i=i+1; gsymm(i) = group_merge(gall([1,3]));
        i=i+1; gsymm(i) = group_merge(gall([2,4]));
        
        i=i+1; gsymm(i) = group_merge(gall([5,7]));
        i=i+1; gsymm(i) = group_merge(gall([6,8]));
        
        i=i+1; gsymm(i) = group_merge(gall([9, 11]));
        i=i+1; gsymm(i) = group_merge(gall([10,12]));
        
        mygroup = gsymm;
    end
    
    %
    
    % Add xdata, legends, xlims, etc. - Need to do this after running
    % group_merge!
    for i = 1:length(mygroup);
        mygroup(i).xdata = out.group(1).xdata;
        mygroup(i).xlims_desired = [0 120];
    end
    
    
    for i = 1:length(mygroup); [mygroup(i).cells, mygroup(i).numcells]= get_grouped_cells(mygroup(i),[sp]); end

    [mygroup] = remove_bads_cells(out.bad_any,mygroup);                 % Remove bad cells from cell lists

    for i = 1:length(mygroup)
        [mygroup(i).data] = get_grouped_data(mygroup(i).ctgs,mygroup(i).cells,pls_sort);
        
        [mygroup(i).metadata.mypairs] = get_grouped_data(ones(size(mygroup(i).ctgs)),mygroup(i).cells,out.mypairs')';
        
        [mygroup(i).datastats, mygroup(i).freqband_stats] = calc_pls_stats(out.abscissa,mygroup(i).data,freqband_stats,'do_mean_ctgs',1);
    end
    

    %mygroup = grouppairs_merge(mygroup,perm2pls_dophi);

    %i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup(:)],[],opts_PM3Dcs,[]);
    %
    if plot_on
        % Estimate ylims
        temp = [];
        mudatastats = arrayfun(@(s) mean(s.datastats),mygroup);
        mylims = [min(mudatastats)*0.9 max(mudatastats)*1.1];
        clear temp
        
        % Xlims
        myxlims = [0 8];
        
        % Do plot
        
        if do_AB_symmetry
%             i=0; i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup(1:2)],[],opts_PM3Dcs,[]); title(mytitles{i})
%             i=0; i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup(3:4)],[],opts_PM3Dcs,[]); title(mytitles{i})
%             i=0; i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup(5:6)],[],opts_PM3Dcs,[]); title(mytitles{i})
            
            opts_PSC.hmask = false(6);
%             i=i+1; figsm; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(:)],opts_PSC);
            
            mygroupp = grouppairs_merge(mygroup,perm2pls_dophi);
            
            for i = 1:length(mygroupp)
                mygroupp(i).legend = ['' mytitles{i}];
            end
            
            opts_PSC.hmask = true(3);
            i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroupp(:)],opts_PSC);
            
            
        end
    end
    
%     i=i+1; figsm; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(1:7)],opts_PSC); title('Sch A Rel'); ylim(mylims); title(mytitles{i});
%     i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(8:14)],opts_PSC); title('Sch B Rel');% ylim([0.45 0.65])
%     i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(15:21)],opts_PSC); title('Sch A Irrel');% ylim([0.45 0.65])
%     i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(22:28)],opts_PSC); title('Sch B Irrel'); %ylim([0.45 0.65])
    % 
    % groupls = out.group;
    % plsls = out.pls_stats;
    % bad_ls = out.bad_any;

    %
    % clear out
end