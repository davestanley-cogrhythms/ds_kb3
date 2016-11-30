function [wrkspc_buffer, group] = Fig_8c_Bndry_vs_Ctgs2_scatter(wrkspc_buffer,curr_stage_sfc,freqband_stats,mygr,opts_exclude,opts_pls,opts_perm,do_abs)

    do_units = 1;

    % Load ctgs pls
    sfc_mode = 22.4014111;
    if do_units; sfc_mode = 52.7000010;end
    [wrkspc_buffer, out] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
    vars_pull(out);
    pls_ctg = pls; pls_stats_ctg = pls_stats;
    ba_ctg = out.bad_any;
    clear out

    % Load boundary pls
    sfc_mode = 22.4514111;
    if do_units; sfc_mode = 52.7500010; end
    [wrkspc_buffer, out] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
    vars_pull(out);
    pls = pls * -1; pls_stats = pls_stats * -1;
    ba = out.bad_any;
    clear out
    
%     Ctgsetli mode 5
%             % New, implemented August 12
%     ctgsetli = [ [t100 & schA] [t50 & schA] [t100 & schB] [t50 & schB] ...
%             [t100_irr & schA] [t50_irr & schA] [t100_irr & schB] [t50_irr & schB] ...
%         ];

    
    % Load ctgs sp
    sfc_mode = 22.4014111;
    if do_units; sfc_mode = 52.7000010; end
    [wrkspc_buffer, out] = load_pr(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,bad_any,opts_perm);
    sp = out.sig_cells; clear out

    % Build group
    group = get_grouped_cells_and_data(mygr,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D);
    
    % Setup variables
    
    % Sch A ctgs vs Sch A Bndry
    figure;
    x = pls_stats_ctg(~ba,1); y = pls_stats(~ba,1);
    if do_abs; x=abs(x);y=abs(y); end
    [x2,y2] = reduce_xy_plusminus(x,y,opts_perm.split_plusminus);
    plot(x,y,'k.');
    hold on; plott_fit(x2,y2,'kx'); % Only fit to data points matching plusminus config.
    hold on; plot_scattergroups(x,y,group([1,4]),~ba); xlabel('Ctg1-Ctg2'); ylabel('SchA B-NB');
    
    % Sch A ctgs vs Sch B Bndry
    figure;
    x = pls_stats_ctg(~ba,1); y = pls_stats(~ba,2);
    if do_abs; x=abs(x);y=abs(y); end
    [x2,y2] = reduce_xy_plusminus(x,y,opts_perm.split_plusminus);
    plot(x,y,'k.');
    hold on; plott_fit(x2,y2,'kx'); % Only fit to data points matching plusminus config.
    hold on; plot_scattergroups(x,y,group([1,4]),~ba); xlabel('Ctg1-Ctg2'); ylabel('SchB B-NB');
    
    % Sch B ctgs vs Sch A Bndry
    figure;
    x = pls_stats_ctg(~ba,2); y = pls_stats(~ba,1);
    if do_abs; x=abs(x);y=abs(y); end
    [x2,y2] = reduce_xy_plusminus(x,y,opts_perm.split_plusminus);
    plot(x,y,'k.');
    hold on; plott_fit(x2,y2,'kx'); % Only fit to data points matching plusminus config.
    hold on; plot_scattergroups(x,y,group([1,4]),~ba); xlabel('Ctg3-Ctg4'); ylabel('SchA B-NB');
    
    % Sch B ctgs vs Sch B Bndry
    figure;
    x = pls_stats_ctg(~ba,2); y = pls_stats(~ba,2);
    if do_abs; x=abs(x);y=abs(y); end
    [x2,y2] = reduce_xy_plusminus(x,y,opts_perm.split_plusminus);
    plot(x,y,'k.');
    hold on; plott_fit(x2,y2,'kx'); % Only fit to data points matching plusminus config.
    hold on; plot_scattergroups(x,y,group([1,4]),~ba); xlabel('Ctg3-Ctg4'); ylabel('SchB B-NB');

    
    
end

