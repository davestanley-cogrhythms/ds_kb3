function [wrkspc_buffer, group] = Fig_6e_rel_vs_irrel(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,opts_PM3Dcs,plot_on)

    do_diff = 0;

    % % Load PLS
    [wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
    pls = out_pls.pls;
    pls_stats = out_pls.pls_stats;
    abscissa = out_pls.abscissa;
    bad_any = out_pls.bad_any;
    funames1D = out_pls.funames1D;
    mypairs = out_pls.mypairs;
    group0 = out_pls.group;
    clear out_pls

    %
    % % Rearrange PLS_perm to be relevant minus irrelevant
    group = group0;
    pls = abs(pls); pls_stats = abs(pls_stats);
    pls = pls(:,:,[1,3,2,4,5]);
    pls_stats = pls_stats(:,[1,3,2,4,5]);
    
    %
    i=0;
    i=i+1; group(i).legend = 'C/D Rel.';
    i=i+1; group(i).legend = 'C/D Irrel.';
    i=i+1; group(i).legend = 'F/T Rel.';
    i=i+1; group(i).legend = 'F/T Irrel.';
    i=i+1; group(i).legend = 'Cue';


    % Setup and load group

    sp = ~bad_any(:);
    group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D);

    if do_diff
        iscirc = opts_pls.perm2pls_dophi;
        group = grouppairs_diff(group,iscirc);
    else
        group = group(1:4);
    end

    %
    % Plot
    if plot_on
        figure;
        inds = 1:length(group);
        [h1] = plot_matrix3D_custstruct([],group(inds),[],opts_PM3Dcs);
    end


    
end

