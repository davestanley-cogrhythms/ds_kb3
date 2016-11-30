function [wrkspc_buffer, group] = Fig_13c0_Units_vs_Ctg(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,perm_mode,curr_stage_sp,freqband_stats_perm,opts_perm,plot_on_func,do_units)

    if ~exist('do_units','var'); do_units = 0; end

    % Load pls
    [wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
    pls = out_pls.pls;
    pls_stats = out_pls.pls_stats;
    abscissa = out_pls.abscissa;
    bad_any = out_pls.bad_any;
    funames1D = out_pls.funames1D;
    mypairs = out_pls.mypairs;
    group0 = out_pls.group;
    clear out_pls

    % Load sp's
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm,opts_exclude);
    sp = out_perm.sig_cells;
    [sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md);

    % Create group template
    mycrit = [2*ones(1,size(sp,2))];
    grt = group0(1);
    grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;
    if do_units; grt.xlims_desired = [-1.5 2.2]; end

    
    % Run a simple test
    clear group
    i=0;
    i=i+1; group(i)=grt; group(i).criteria(1:2)=[1 0]; group(i).ctgs=1; % Ctg1-2 deciders
    i=i+1; group(i)=grt; group(i).criteria(1:2)=[1 0]; group(i).ctgs=2; % Ctg1-2 deciders
    i=i+1; group(i)=grt; group(i).criteria(1:2)=[0 0]; group(i).ctgs=[1]; % Ctg1-2 non-deciders
    i=i+1; group(i)=grt; group(i).criteria(1:2)=[0 1]; group(i).ctgs=2;   % Ctg1-2 deciders
    i=i+1; group(i)=grt; group(i).criteria(1:2)=[0 1]; group(i).ctgs=1;   % Ctg1-2 non-deciders
    i=i+1; group(i)=grt; group(i).criteria(1:2)=[0 0]; group(i).ctgs=[2];   % Ctg1-2 non-deciders
    

    % Calculate legend entries
    group = group.query_legend(group0);

    

    %for i = 1:4; group(i).interpreter = 'latex'; end
    
    if do_units;
        %pername = 'Esb.';
        %pername = 'Units';
        pername = 'Ens.';
    else pername = 'Units';
    end
    
    i=0;
    i=i+1;  group(i).legend = 'A Pr';
    i=i+1;  group(i).legend = 'A NP';
    i=i+1;  group(i).legend = 'ND A';
    i=i+1;  group(i).legend = 'B Pr';
    i=i+1;  group(i).legend = 'B NP';
    i=i+1;  group(i).legend = 'ND B';
    
    group1=group_merge(group(1),group(4));
    group2=group_merge(group(2),group(5));
    group3=group_merge(group(3),group(6));
    
    group = [group1 group2 group3];
    
    
    % Load data into groups
    group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D);
    
    
    if plot_on_func
        figure;
        inds = 1:length(group);
        [h1] = plot_matrix3D_custstruct([],group(inds));
    end

    % Test bargraph

    if plot_on_func
        i=i+1; fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group);
    end



end