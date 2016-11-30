

function mygroup = Fig_7_Ctg_Selectivity(wrkspc_buffer,out,perm_mode,perm_mode2,curr_stage_sfc,freqband_stats,perm2pls_dophi)

    plot_on = 1;
    
    % Plot matrix settings
    opts_PM3Dcs.paperfig_mode = 1;
    opts_PM3Dcs.do_sgolay = 0;
    opts_PM3Dcs.stats_mode = 0;

    % Statistics settings
    opts_PPS_temp = Opts_Perm;
    opts_PPS_temp.do_bh0=1;
    opts_PPS_temp.do_phi=perm2pls_dophi;
    opts_PPS_temp.split_plusminus=0;

    
    % Sort into preferred and non-preferred
    Cave1 = wrkspc_buffer.(['per_' mode2modename(perm_mode2)]).(stagename(curr_stage_sfc)).Cave1;
    Cave2 = wrkspc_buffer.(['per_' mode2modename(perm_mode2)]).(stagename(curr_stage_sfc)).Cave2;


    [IA] = ind_pls(Cave1(:,:,1),Cave2(:,:,1),out.abscissa,freqband_stats);  % Sch A
    [IB] = ind_pls(Cave1(:,:,2),Cave2(:,:,2),out.abscissa,freqband_stats);  % Sch B

    sz = size(out.pls);
    pls_sort = zeros([sz(1),sz(2),8]);
    pls_sort(:,:,[1:2]) = sort_index_pls(out.pls(:,:,1),out.pls(:,:,2),IA);      % 60% A
    pls_sort(:,:,[3:4]) = sort_index_pls(out.pls(:,:,3),out.pls(:,:,4),IB);      % 60% B
    pls_sort(:,:,[5:6]) = sort_index_pls(out.pls(:,:,9),out.pls(:,:,10),IA);     % 100% A
    pls_sort(:,:,[7:8]) = sort_index_pls(out.pls(:,:,11),out.pls(:,:,12),IB);    % 100% B

    mylegend = {'SchA 60','SchB 60','SchA 100','SchB 100'};

    %

    % Get statistics
    s = wrkspc_buffer.(['per_' mode2modename(perm_mode)]).(stagename(curr_stage_sfc));
    [sp] = plot_permute_scatter_fast(s,out.bad_any,freqband_stats,opts_PPS_temp);
    sp = sp(:,[1,2,5,6]);
    clear s

    %
    % Build group
    myN_criteria = 4;
    gr_template = Grp;
    gr_template.criteria = [2*ones(1,myN_criteria)]; gr_template.criteria_alt = []; gr_template.criteria_sfc = []; gr_template.ctgs = 1;
    gr_template.xlims_desired = [0 120]; gr_template.xdata = out.group(1).xdata;


    % indA = [1 0 1 0];
    % indB = [0 1 0 1];
    % indA = [1 0 2 2];
    % indB = [0 1 2 2];
    % indA = [2 2 1 0];
    % indB = [2 2 0 1];


    indA = [2 2 2 2]; indB = indA;

    clear group
    i=0;
    i=i+1; mygroup(i) = gr_template; mygroup(i).criteria=[indA]; mygroup(i).ctgs = [1]; mygroup(i).legend = 'Sch A Pref 60';
    i=i+1; mygroup(i) = gr_template; mygroup(i).criteria=[indA]; mygroup(i).ctgs = [2]; mygroup(i).legend = 'Sch A Non-pref 60';
    i=i+1; mygroup(i) = gr_template; mygroup(i).criteria=[indA]; mygroup(i).ctgs = [5]; mygroup(i).legend = 'Sch A Pref 100';
    i=i+1; mygroup(i) = gr_template; mygroup(i).criteria=[indA]; mygroup(i).ctgs = [6]; mygroup(i).legend = 'Sch A Non-pref 100';
    i=i+1; mygroup(i) = gr_template; mygroup(i).criteria=[indB]; mygroup(i).ctgs = [3]; mygroup(i).legend = 'Sch B Pref 60';
    i=i+1; mygroup(i) = gr_template; mygroup(i).criteria=[indB]; mygroup(i).ctgs = [4]; mygroup(i).legend = 'Sch B Non-pref 60';
    i=i+1; mygroup(i) = gr_template; mygroup(i).criteria=[indB]; mygroup(i).ctgs = [7]; mygroup(i).legend = 'Sch B Pref 100';
    i=i+1; mygroup(i) = gr_template; mygroup(i).criteria=[indB]; mygroup(i).ctgs = [8]; mygroup(i).legend = 'Sch B Non-pref 100';

    for i = 1:length(mygroup); [mygroup(i).cells, mygroup(i).numcells]= get_grouped_cells(mygroup(i),[sp]); end

    [mygroup] = remove_bads_cells(out.bad_any,mygroup);                 % Remove bad cells from cell lists

    for i = 1:length(mygroup) [mygroup(i).data] = get_grouped_data(mygroup(i).ctgs,mygroup(i).cells,pls_sort); end

    mygroup = grouppairs_merge(mygroup,perm2pls_dophi);

    if plot_on
        i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup(1:end/2)],[],opts_PM3Dcs,[]);
        i=i+1; figure; [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[mygroup(end/2+1:end)],[],opts_PM3Dcs,[]);
    end

    % 
    % groupls = out.group;
    % plsls = out.pls_stats;
    % bad_ls = out.bad_any;

    %
    clear out
end