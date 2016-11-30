function mygroup = Fig_9_Royfig(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)
    %% Setup
    % Settings
    plot_on = 1;
    do_AB_symmetry = 1;
    
    cell_group_mode = 2;        %1-Do all cells; 2-Only ctg deciders; 3-only boundary deciders
        if cell_group_mode == 2 || cell_group_mode == 3; get_statistics = 1;
        else get_statistics = 0;
        end
    
    %
    % Plot matrix settings
    opts_PM3Dcs.paperfig_mode = 1;
    opts_PM3Dcs.do_sgolay = 0;
    opts_PM3Dcs.stats_mode = 0;

    % Plot bar settings
    opts_PSC.paperfig_mode = 1;
    opts_PSC.hmask = false(7);
    
    % Statistics settings
    opts_PPS_temp.do_bh0=1;
    opts_PPS_temp.do_phi=perm2pls_dophi;
    opts_PPS_temp.split_plusminus = 0;
    
    
    

    %% Sorting into pref / non-pref
    % Sort into preferred and non-preferred
    Cave1 = wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)).Cave1;
    Cave2 = wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc)).Cave2;


    [IA] = ind_pls(Cave1(:,:,1),Cave2(:,:,1),out.abscissa,freqband_stats);  % Sch A
    [IB] = ind_pls(Cave1(:,:,2),Cave2(:,:,2),out.abscissa,freqband_stats);  % Sch B
    warning('Maybe redo for irrelevant deciders.');
    

    sz = size(out.pls);
    pls_sort = zeros([sz(1),sz(2),28]);
    
    del_p=8;
    del_ps = 7;
    % Sch A Relevant [1,2]; arrange into "Roy" type structure.
    pls_sort(:,:,[3+0*del_ps,5+0*del_ps]) = sort_index_pls(out.pls(:,:,1+0*del_p),out.pls(:,:,2+0*del_p),IA);  % 60%
    pls_sort(:,:,[2+0*del_ps,6+0*del_ps]) = sort_index_pls(out.pls(:,:,1+1*del_p),out.pls(:,:,2+1*del_p),IA);  % 80%
    pls_sort(:,:,[1+0*del_ps,7+0*del_ps]) = sort_index_pls(out.pls(:,:,1+2*del_p),out.pls(:,:,2+2*del_p),IA);  % 100%
    pls_sort(:,:,[4+0*del_ps]) = out.pls(:,:,del_p*3+1);
    
    % Sch B Relevant [3,4]
    pls_sort(:,:,[3+1*del_ps,5+1*del_ps]) = sort_index_pls(out.pls(:,:,3+0*del_p),out.pls(:,:,4+0*del_p),IB);  % 60%
    pls_sort(:,:,[2+1*del_ps,6+1*del_ps]) = sort_index_pls(out.pls(:,:,3+1*del_p),out.pls(:,:,4+1*del_p),IB);  % 80%
    pls_sort(:,:,[1+1*del_ps,7+1*del_ps]) = sort_index_pls(out.pls(:,:,3+2*del_p),out.pls(:,:,4+2*del_p),IB);  % 100%
    pls_sort(:,:,[4+1*del_ps]) = out.pls(:,:,del_p*3+2);
    
    % Sch A IrrRelevant (i.e. cued for doing fat/thin, but pls is measuring
    % cat/dog
    pls_sort(:,:,[3+2*del_ps,5+2*del_ps]) = sort_index_pls(out.pls(:,:,5+0*del_p),out.pls(:,:,6+0*del_p),IA);  % 60%
    pls_sort(:,:,[2+2*del_ps,6+2*del_ps]) = sort_index_pls(out.pls(:,:,5+1*del_p),out.pls(:,:,6+1*del_p),IA);  % 80%
    pls_sort(:,:,[1+2*del_ps,7+2*del_ps]) = sort_index_pls(out.pls(:,:,5+2*del_p),out.pls(:,:,6+2*del_p),IA);  % 100%
    pls_sort(:,:,[4+2*del_ps]) = out.pls(:,:,del_p*3+4);        % Sch B cue + 50% irrelevant

    % Sch B IrrRelevant (i.e. cued for doing cat/dog, but pls is measuring
    % fat/thin
    pls_sort(:,:,[3+3*del_ps,5+3*del_ps]) = sort_index_pls(out.pls(:,:,7+0*del_p),out.pls(:,:,8+0*del_p),IB);  % 60%
    pls_sort(:,:,[2+3*del_ps,6+3*del_ps]) = sort_index_pls(out.pls(:,:,7+1*del_p),out.pls(:,:,8+1*del_p),IB);  % 80%
    pls_sort(:,:,[1+3*del_ps,7+3*del_ps]) = sort_index_pls(out.pls(:,:,7+2*del_p),out.pls(:,:,8+2*del_p),IB);  % 100%
    pls_sort(:,:,[4+3*del_ps]) = out.pls(:,:,del_p*3+3);        % Sch A cue + 50% irrelevant

    mylegend = {'100','80','60','50','40','20','0'};

    %% Get statistics
    if get_statistics
        % Do sp for perm_mode (ctg sensitivity)
        s = wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc));
        [sp1] = plot_permute_scatter_fast(s,out.bad_any,freqband_stats,opts_PPS_temp);
        clear s
        
        % Do sp for perm_mode2 (bndry sensitivity)
        s = wrkspc_buffer.(['per_' mode2modename(perm_mode2)]).(stagename(curr_stage_sfc));
        [sp2] = plot_permute_scatter_fast(s,out.bad_any,freqband_stats,opts_PPS_temp);
        clear s
        sp = [sp1 sp2];
    else
        Ncrit = 9;
        sp = false(size(pls_sort,2),Ncrit);
    end

    %% Setup groups
    % Build group template
    mycrit = [2*ones(1,size(sp,2))];
    gr_template = Grp;
    gr_template.criteria = mycrit; gr_template.criteria_alt = []; gr_template.criteria_sfc = []; gr_template.ctgs = 1;
    gr_template.xlims_desired = [0 120]; gr_template.xdata = out.group(1).xdata;

    clear gall
    switch cell_group_mode
        case 1  % Anything
            mycrit = [2,2,2,2,2,   2,2,2,2];
            for i = 1:size(pls_sort,3)
                gall(i) = gr_template; gall(i).criteria=[mycrit]; gall(i).ctgs = i; 
            end
        case 2  % Deciders
            mycritA = [1,0,2,2,2,   2,2,2,2];
            mycritB = [0,1,2,2,2,   2,2,2,2];
            
            % Sch A preferred
            for i = 1:7; gall(i) = gr_template; gall(i).criteria=[mycritA]; gall(i).ctgs = i; end
            
            % Sch B preferred
            for i = 8:14; gall(i) = gr_template; gall(i).criteria=[mycritB]; gall(i).ctgs = i; end

            % Sch A Irr (i.e. Scheme A categorizations when Scheme B is cued - i.e. scheme A irrelevant)
            for i = 15:21; gall(i) = gr_template; gall(i).criteria=[mycritA]; gall(i).ctgs = i; end
            
            % Sch B Irr (i.e. Scheme B categorizations when Scheme A is cued - i.e. scheme B irrelevant)
            for i = 22:28; gall(i) = gr_template; gall(i).criteria=[mycritB]; gall(i).ctgs = i; end
            
            % Non-preferreds
            gall_np = gall(1:14);
            for i = 1:7; gall_np(i).criteria = mycritB; % Sch B NP
            end
            
            for i = 8:14; gall_np(i).criteria = mycritA; % Sch A NP
            end
            gall_np = gall_np([8:14,1:7]);  % Swap so Sch A non-preferred comes first
            
            % Non-preferrds become 29-35 (SchA) and 36-42 (SchB)
            gall = [gall gall_np];
            
            
        case 3
                % Boundary electrodes
            mycritA = [2,2,2,2,2,   1,0,2,2];
            mycritB = [2,2,2,2,2,   0,1,2,2];
            mycrit = [repmat(mycritA,7,1); repmat(mycritB,7,1); repmat(mycritA,7,1); repmat(mycritB,7,1); ];
            
            for i = 1:size(mycrit,1);
                gall(i) = gr_template; gall(i).criteria=[mycrit(i,:)]; gall(i).ctgs = i;
            end
    end
    
    % Build legend & titles
    for i = 1:length(gall); gall(i).legend = mylegend{mod(i-1,length(mylegend))+1}; end
    mytitles={'A Preferred Rel','B Preferred Rel','A Pref Irrel','B Pref Irrel','A Non-pref Rel','B Non-pref Rel'};
    mygroup = gall;
    
    
    if do_AB_symmetry
        clear gsymm
        for i = 1:size(pls_sort,3)/4
            gsymm(i) = group_merge(gall([i,i+del_ps]));          % Merge 1..7 and 8..14  - Preferred Relevant
            gsymm(i+7) = group_merge(gall([i+2*del_ps,i+3*del_ps])); % Merge 15...21 and 22...28 - Preferred irrelevant
            if cell_group_mode == 2 % Non-preferred relevant case is only when we have cells with a preference!
                gsymm(i+14) = group_merge(gall([i+4*del_ps,i+5*del_ps])); % Merge 29...35 and 37...42 - Non-preferred relevant (if applicable)
            end
        end
        %Build legend & titles
        for i = 1:length(gsymm); gsymm(i).legend = mylegend{mod(i-1,length(mylegend))+1}; end
        mytitles={'Preferred Rel','Pref Irrel','Non-pref Rel'};
        mygroup = gsymm;
        
    end
    
    %
    for i = 1:length(mygroup); [mygroup(i).cells, mygroup(i).numcells]= get_grouped_cells(mygroup(i),[sp]); end

    [mygroup] = remove_bads_cells(out.bad_any,mygroup);                 % Remove bad cells from cell lists

    for i = 1:length(mygroup)
        [mygroup(i).data] = get_grouped_data(mygroup(i).ctgs,mygroup(i).cells,pls_sort);
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
        
        opts_PSC.hmask = false(7);
        if do_AB_symmetry
            i=0;
            i=i+1; figsm; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(1:7)],opts_PSC); ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(1).numcells))]);  xlim(myxlims);
            i=i+1; figsm; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(8:14)],opts_PSC); ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(8).numcells))]); xlim(myxlims);
            if length(mygroup) > 14     % We're doing ctg deciders (cell group mode == 2)
                i=i+1; figsm; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(15:21)],opts_PSC); ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(15).numcells))]); xlim(myxlims);
            end
        else
            
            
            i=0;
            i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(1:7)],opts_PSC);ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(1).numcells))]);  xlim(myxlims);
            i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(8:14)],opts_PSC); ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(8).numcells))]);  xlim(myxlims);
            i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(15:21)],opts_PSC); ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(15).numcells))]);  xlim(myxlims);
            i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct([mygroup(22:28)],opts_PSC); ylim(mylims); title(mytitles{i}); legend(['N = ' num2str(sum(mygroup(22).numcells))]);  xlim(myxlims);
        end
    end
    
    % 
    % groupls = out.group;
    % plsls = out.pls_stats;
    % bad_ls = out.bad_any;

    %
    % clear out
end