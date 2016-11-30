function [wrkspc_buffer, group, group_m ] = Fig_4c2_Bndry_redo(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,opts_PM3Dcs,opts_PM3Dsp,plot_on_func)

    grouppairs_merge_operation = 1;
    
    % Load pls
    [wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
    pls = out_pls.pls;
    pls_stats = out_pls.pls_stats;
    abscissa = out_pls.abscissa;
    abscissa2 = out_pls.abscissa2;
    bad_any = out_pls.bad_any;
    funames1D = out_pls.funames1D;
    mypairs = out_pls.mypairs;
    group0 = out_pls.group;
    clear out_pls

    % Build group
    % Load sp's
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm, opts_exclude);
    sp = out_perm.sig_cells;

    % Load bads perm
    [bad_any_perm] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);

    % Map sp's as needed
    [sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any);
    
    group = group0;
    sp = ~bad_any(:);

    % Load data into groups
    group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);
    data_type = mode2datatype(sfc_mode);
    for i = 1:length(group); group(i).data_name = data_type; end
    
    group_m = grouppairs_merge(group([2,1,4,3,6,5,8,7,10,9,12,11]),opts_pls.perm2pls_dophi,grouppairs_merge_operation);        % Swap 100 and 50%'s
    if grouppairs_merge_operation
        for i = 1:length(group_m); group_m(i).data_name = '% change'; end
    else
        for i = 1:length(group_m); group_m(i).data_name = 'Diff'; end
    end

    % Test Plot groups, at last

    if plot_on_func
        opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};


        inds = 1:length(group);
        figure;
        [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);

    %     figure;
    %     inds = 1:4;
    %     [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D);
    %     figure;
    %     inds = 5:8;
    %     [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D);
    end

    % Test bargraph

    if plot_on_func
        fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group);
    end

end
