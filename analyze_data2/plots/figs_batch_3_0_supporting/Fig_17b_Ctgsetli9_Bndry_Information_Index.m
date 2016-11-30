function [wrkspc_buffer, gall] = Fig_17b_Ctgsetli9_Bndry_Information_Index(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode_sort,perm_mode_sp,curr_stage_sp,freqband_stats_perm,opts_perm,do_unitsCI_vs_electBNDRY,plot_on_func)


    %%
    sort_on = 1; % Setting to 0 completely disables all sorting
    do_sort_irr = 1; % 1-process irrelevants separately; 0-sort everything based on relevant category(maybe biased)

    [wrkspc_buffer, out_pls, out_sort, out_perm, sp] = Fig_9_general_load_and_arrange_CI(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode_sort,perm_mode_sp,curr_stage_sp,freqband_stats_perm,opts_perm,do_unitsCI_vs_electBNDRY);
    pls = out_pls.pls;
    pls_stats = out_pls.pls_stats;
    abscissa = out_pls.abscissa;
    bad_any = out_pls.bad_any;
    funames1D = out_pls.funames1D;
    mypairs = out_pls.mypairs;
    group0 = out_pls.group;
    
    sp = zeros(size(pls,2), size(sp,2));  % % % %% %% % % %% % OVERWWRITE SP!!

    if ~opts_pls.collapse_pls_to_days
        pls_sort = get_pls_sort(out_pls, out_sort, sort_on, do_sort_irr, freqband_stats);
    else
        pls_sort = pls;
    end

    %%

    del_ps = 7;
    lr = 1:3;   % Left categories
    rr = 5:7;   % Right categories
    pls_ctgs = cat(3,  mean(pls_sort(:,:,lr + 0*del_ps),3) - mean(pls_sort(:,:,rr + 0*del_ps),3), ...
                       mean(pls_sort(:,:,lr + 1*del_ps),3) - mean(pls_sort(:,:,rr + 1*del_ps),3), ...
                       mean(pls_sort(:,:,lr + 2*del_ps),3) - mean(pls_sort(:,:,rr + 2*del_ps),3), ...
                       mean(pls_sort(:,:,lr + 3*del_ps),3) - mean(pls_sort(:,:,rr + 3*del_ps),3)  ...
                       );

    br = 4;
    nbr = [1,7];
    pls_bnds = cat(3,  mean(pls_sort(:,:,br + 0*del_ps),3) - mean(pls_sort(:,:,nbr + 0*del_ps),3), ...
                       mean(pls_sort(:,:,br + 1*del_ps),3) - mean(pls_sort(:,:,nbr + 1*del_ps),3), ...
                       mean(pls_sort(:,:,br + 2*del_ps),3) - mean(pls_sort(:,:,nbr + 2*del_ps),3), ...
                       mean(pls_sort(:,:,br + 3*del_ps),3) - mean(pls_sort(:,:,nbr + 3*del_ps),3)  ...
                       );
    pls_ii = ((pls_bnds) - abs(pls_ctgs)) ./ (abs(pls_bnds) + abs(pls_ctgs));



    % Group stuff
    mycrit = [2*ones(1,size(sp,2))];
    grt = group0(1);
    grt.criteria = []; grt.criteria_alt = []; grt.criteria_sfc = [];
    grt.xdata = group0(1).xdata;
    grt.data_name = ['Category Index'];
    grt.interpreter = 'latex';

    clear gall
    for i = 1:size(pls_ii,3)
        gall(i)=grt; gall(i).criteria=[mycrit]; gall(i).ctgs = i;
    end


    i=0;
    i=i+1; gall(i).legend='BII Cat/Dog';
    i=i+1; gall(i).legend='BII Goc/Tad';
    i=i+1; gall(i).legend='BII Cat/Dog irr';
    i=i+1; gall(i).legend='BII Goc/Tad irr';

    
    gall = get_grouped_cells_and_data(gall,sp,pls_ii,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D);

    if plot_on_func
        figure;
        [h1] = plot_matrix3D_custstruct([],gall);
    end





end