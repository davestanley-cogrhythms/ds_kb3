function [wrkspc_buffer, group, pls, pls_stats, bad_any] = load_pls_sp_groupdata(wrkspc_buffer,sfc_mode,perm_mode,stage_sfc,stage_perm,freqband_stats,mygr,opts_exclude,opts_pls,opts_perm)

    % Load boundary pls
    %sfc_mode = 22.4514111;
    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,stage_sfc,freqband_stats,opts_exclude,opts_pls);
    %vars_pull(out_pls);
    bad_any = out_pls.bad_any;
    pls = out_pls.pls;
    pls_stats = out_pls.pls_stats;
    abscissa = out_pls.abscissa;
    mypairs = out_pls.mypairs;
    funames1D = out_pls.funames1D;
    
    if mode_subgroups(2) == 5 && ~opts_pls.perm2pls
        pls = pls * -1; pls_stats = pls_stats * -1;
    end
    
    

    % Load ctgs sp
    %perm_mode = 22.4014111;
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,stage_perm,freqband_stats,bad_any,opts_perm);
    sp = out_perm.sig_cells; clear out_perm

    % Build group
    group = get_grouped_cells_and_data(mygr,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D);
    
    
end