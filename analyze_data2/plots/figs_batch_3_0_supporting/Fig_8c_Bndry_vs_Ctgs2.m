function [wrkspc_buffer, group] = Fig_8c_Bndry_vs_Ctgs2(wrkspc_buffer,curr_stage_sfc,freqband_stats,mygr,opts_exclude,opts_pls,opts_perm)

    do_units = 1;

    % Load boundary pls
    sfc_mode = 22.4514111;
    if do_units; sfc_mode = 52.7500010;end
    [wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
    vars_pull(out_pls);
    pls = pls * -1; pls_stats = pls_stats * -1;
    clear out_pls


    % Load ctgs sp
    sfc_mode = 22.4014111;
    if do_units; sfc_mode = 52.7000010; end
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,bad_any,opts_perm);
    sp = out_perm.sig_cells; clear out_perm

    % Build group
    group = get_grouped_cells_and_data(mygr,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D);
    
    
end