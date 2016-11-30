function [wrkspc_buffer, out_pls, out_sort, out_perm, sp] = Fig_9_general_load_and_arrange_CI(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode_sort,perm_mode_sp,curr_stage_sp,freqband_stats_perm,opts_perm,do_unitsCI_vs_electBNDRY)

    % Setup unitsCI_vs_electBNDRY parameters
    opts_perm.chosen_quantile = 0.2;
    opts_perm.upper_quantile = 0;       % Do lower quantile, since boundary trials are 100% - 50%

    % % Load PLS
    [wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
    pls = out_pls.pls;
    pls_stats = out_pls.pls_stats;
    abscissa = out_pls.abscissa;
    bad_any = out_pls.bad_any;
    funames1D = out_pls.funames1D;
    mypairs = out_pls.mypairs;
    group0 = out_pls.group;

    % % Load perm mode for sort
    [wrkspc_buffer, out_sort] = load_pr(wrkspc_buffer,perm_mode_sort,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm,opts_exclude);
    
    % % Load perm mode for sp
    if do_unitsCI_vs_electBNDRY;
        opts_perm.do_quantiles_mode = 1;
    end
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode_sp,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm,opts_exclude);
    sp = out_perm.sig_cells;

    % Load bads perm
    [bad_any_perm] = load_bads(perm_mode_sp,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);

    % Map sp's as needed
    [sp] = map_sp(perm_mode_sp, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any);

end


