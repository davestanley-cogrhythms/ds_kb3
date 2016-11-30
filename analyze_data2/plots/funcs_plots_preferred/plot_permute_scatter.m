
function [out] = plot_permute_scatter(s_perm,s_boot, bad_any,freqband_stats, opts)
    %% plot_SchAvsB_scatter - Plot significant differences AMPLITUDE across SchA and SchB
    permmode = 1;               % Permmode: 1-Significant difference between Schs; 2-no significant difference between schemes
    
    if nargin < 5
        opts = [];
    end
    
    % Plot switches
    if ~isfield(opts,'plot_Sch_sens_spect'); opts.plot_Sch_sens_spect = false; end
    if ~isfield(opts,'plot_Sch_sens_scatter'); opts.plot_Sch_sens_scatter = true; end
    if ~isfield(opts,'plot_Sch_sens_buschman'); opts.plot_Sch_sens_buschman = false; end
    if ~isfield(opts,'plot_Sch_sens_cumulative'); opts.plot_Sch_sens_cumulative = false; end
    if ~isfield(opts,'hide_all_plots'); opts.hide_all_plots = 0; end
    if ~isfield(opts,'hide_boot_results'); opts.hide_boot_results = 0; end
    if ~isfield(opts,'do_bh0'); opts.do_bh0 = 0; end
    if ~isfield(opts,'do_phi'); opts.do_phi = 0; end
    if ~isfield(opts,'paperfig_mode'); opts.paperfig_mode = 0; end      % Not used so far
        
    % Naming 
    if ~isfield(opts,'sens_ctg1_name'); opts.sens_ctg1_name = 'Sch A'; end
    if ~isfield(opts,'sens_ctg2_name'); opts.sens_ctg2_name = 'Sch B'; end
    
    % Other plotting switches
    if ~isfield(opts,'sens_chosen_condition'); opts.sens_chosen_condition = 1; end
    if ~isfield(opts,'pd0'); opts.pd0 = 30; end
    
    % Pull in options
    vars_pull(opts);
    
    % Apply hide all plots
    plot_Sch_sens_spect = plot_Sch_sens_spect && ~hide_all_plots;
    plot_Sch_sens_scatter = plot_Sch_sens_scatter && ~hide_all_plots;
    plot_Sch_sens_buschman = plot_Sch_sens_buschman && ~hide_all_plots;
    plot_Sch_sens_cumulative = plot_Sch_sens_cumulative && ~hide_all_plots;
    
    chosen_cells = ~bad_any;
    
    percent_change_mode = 1;    % 1-(Preferred - Nonpref) / Nonpref ; 2-(Pref - Nonpref) / Pref
    alpha = 0.05;
    if do_bh0; alpha = 0.20; end
    
    % Get permmode 1 data for adding to graph below
    do_bh = do_bh0;
    s = s_perm;
    if ~do_phi; pvals_Cave_perm = s.pvals_Cave(:,:,sens_chosen_condition);
    else pvals_Cave_perm = s.pvals_phi(:,:,sens_chosen_condition); end
    f = s.f;
    [q4_perm, pvals1, p0] = calc_significant_permutation(f,pvals_Cave_perm,alpha,bad_any,mean(freqband_stats),do_bh);
    q4_perm = q4_perm & chosen_cells;
    sum(q4_perm)
    
    % Code to group significant differences as positive or negative
    ind = find_closest(f,mean(freqband_stats));
    signs = sign(s.Cave1(ind,:,sens_chosen_condition) -s.Cave2(ind,:,sens_chosen_condition));
    q4_permp = q4_perm(:) & signs(:) >= 0;
    q4_permm = q4_perm(:) & signs(:) < 0;
    
    if permmode == 1
        do_bh = do_bh0;
        s = s_perm;
        if ~do_phi; pvals_Cave = s.pvals_Cave(:,:,sens_chosen_condition);
        else pvals_Cave = s.pvals_phi(:,:,sens_chosen_condition); end
        pd = pd0;

    elseif permmode == 2
        do_bh = 1;
        s = s_boot;

        vars_pull(s);

        ksp = horzcat(ks.p);
        ksstat = horzcat(ks.ksstat);
        ztp = horzcat(zerotest.p);
        ztpd = horzcat(zerotest.pd);
        sig = dCave_bs_sig;
        
        
        
        ind = find_closest(f,mean(freqband_stats));
        pd_range = 1:1:500;
        pvals_Cave_range = [];
        
        x = abs(Cave1 - Cave2);
        x = x(ind,:);
        if percent_change_mode==1
            mu = min(Cave1,Cave2);
        elseif percent_change_mode==2
            mu = max(Cave1,Cave2);
        end
        mu = mu(ind,:);
        
        if plot_Sch_sens_cumulative
            for i = 1:length(pd_range)

                pd = pd_range(i);    % percent

                % Implement method to get p value
                delta = mu * (pd/100);                % Reject null hypothesis for changes in Cave less than this value
                %delta = 0.02

                sz = size(x);
                pvals_Cave = 1 - normcdf(delta - x,zeros(sz),sig(ind,:));
                pvals_Cave_range = cat(3,pvals_Cave_range,pvals_Cave);

            end
        end
        
        pd = pd0;    % percent
        x = abs(Cave1 - Cave2);

        % Implement method to get p value
        if percent_change_mode==1
            mu = min(Cave1,Cave2);
        elseif percent_change_mode==2
            mu = max(Cave1,Cave2);
        end
        delta = mu * (pd/100);                % Reject null hypothesis for changes in Cave less than this value
        %delta = 0.02

        sz = size(x);
        pvals_Cave = 1 - normcdf(delta - x,zeros(sz),sig);
        h = pvals_Cave < 0.05;

        %figure; imagesc(h);
        ind = find_closest(f,mean(freqband_stats));
        mysum = sum(h(ind,:));

    end
    

    if ~do_phi
        Cave1 = s.Cave1(:,:,sens_chosen_condition);
        Cave2 = s.Cave2(:,:,sens_chosen_condition);
    else
        Cave1 = s.phi1(:,:,sens_chosen_condition);
        Cave2 = s.phi2(:,:,sens_chosen_condition);
        
        Cave1_orig = s.Cave1(:,:,sens_chosen_condition);    % Always Cave
        Cave2_orig = s.Cave2(:,:,sens_chosen_condition);
    end
    f = s.f;


    
    %[q4_sfc, pvals1] = calc_significant_permutation_old_withplots(s.f,s.dCave,s.Cavep1,s.Cavep2,alpha,bad_any,mean(freqband_stats),do_bh,Cave); sum(q4_sfc)
    [q4_boot, pvals1, p0] = calc_significant_permutation(f,pvals_Cave,alpha,bad_any,mean(freqband_stats),do_bh);
    q4_boot = q4_boot & chosen_cells;
    sum(q4_boot)
    
    % q4_sfc2 = x(ind,:) < delta(ind,:);  % Threshold difference instead of doing stats

    freqrange_str = [num2str(freqband_stats(1)) ' - ' num2str(freqband_stats(2)) ' Hz'];
    var_str = 'FFC';
    Nall_str = ['N=' num2str(sum(chosen_cells))];
    Nnull_str = ['N=' num2str(sum(q4_boot))];
    Nperm_str = ['N=' num2str(sum(q4_perm))];
    
    
    
    
    
    I = ind_pls(Cave1,Cave2,f,freqband_stats);
    Cavepref = sort_index_pls(Cave1,Cave2,I);
    if plot_Sch_sens_spect
        figure; plott_matrix3D(Cavepref(:,q4_boot,:));xlim([0 100]); xlabel('Freq (Hz)'); ylabel(var_str); legend(['Preferred ' Nnull_str],['Non-preferred ' Nnull_str]);
        title(['Preferred vs Non-Preferred']);
    end
    ind = find_closest(f,mean(freqband_stats));
    %Cavef = Cavepref(ind,:,:);
    Cavef = cat(3,Cave1(ind,:,:),Cave2(ind,:,:));
    
    
    % If we're in do_phi mode, recalculate Cavepreferred also using Cave data.
    if do_phi
        I = ind_pls(Cave1_orig,Cave1_orig,f,freqband_stats);
        Cavepreferred_orig = sort_index_pls(Cave1_orig,Cave1_orig,I);
    else
        Cavepreferred_orig = Cavepref(:,:,1);
    end
    

    % Scatterplot 
    if plot_Sch_sens_scatter
        Cavef = squeeze(Cavef);
        figure;
        plot(Cavef(chosen_cells,1),Cavef(chosen_cells,2),'k.');
        if ~hide_boot_results;
            hold on; plot(Cavef(q4_boot,1),Cavef(q4_boot,2),'b.','MarkerSize',20,'LineWidth',2);
        end
        hold on; plot(Cavef(q4_perm,1),Cavef(q4_perm,2),'r.','MarkerSize',20,'LineWidth',2);
        xlabel([var_str ' Sch A (' freqrange_str ')']); ylabel([var_str ' Sch B (' freqrange_str ')']); % title(['Freqrange ' freqrange_str]);
        
        if ~hide_boot_results
            legend(['All data ' Nall_str],['d' var_str ' < ' num2str(pd) '% ' Nnull_str],['d' var_str ' > 0% ' Nperm_str]);
        else
            legend(['All data ' Nall_str],['Sig ' var_str ' > 0% ' Nperm_str]);
        end
    end
    
    % Buschman Plot
    if plot_Sch_sens_buschman
        figure; hold on;
        Cave_preferred = Cavepref(:,:,1); % dim3 index1 is preferred, index2 is non-preferred
        Cave_nonpreferred = Cavepref(:,:,2); % dim3 index1 is preferred, index2 is non-preferred
        dCave = Cave1 - Cave2;
        if ~do_phi
            if percent_change_mode == 1; Cave_percent_change = (dCave) ./ Cave_nonpreferred * 100;
            elseif percent_change_mode == 2; Cave_percent_change = (dCave) ./ Cave_preferred * 100; end
        else
             dCave = get_circ_diff(dCave);  % If in phi mode, only care about difference in phases (normalization is pointless).
             Cave_percent_change = dCave;
        end
        plot(Cavepreferred_orig(ind,chosen_cells),Cave_percent_change(ind,chosen_cells),'k.'); 
        if ~hide_boot_results
            plot(Cavepreferred_orig(ind,q4_boot),Cave_percent_change(ind,q4_boot),'b.','MarkerSize',20);
        end
        plot(Cavepreferred_orig(ind,q4_perm),Cave_percent_change(ind,q4_perm),'r.','MarkerSize',20,'LineWidth',2);
        minmax = [min(Cavepreferred_orig(ind,chosen_cells)) max(Cavepreferred_orig(ind,chosen_cells))];
        plot([minmax],[pd pd],'k:','LineWidth',1)
        plot([minmax],-1*[pd pd],'k:','LineWidth',1)
        xlabel([var_str ' Preferred']); ylabel(['<-- Pref ' sens_ctg2_name ' - d' var_str ' (%) - Pref ' sens_ctg1_name ' -->'])
        title(['Freqrange (' freqrange_str ')']); 
        if ~hide_boot_results
            legend(['All data ' Nall_str],['d' var_str ' < ' num2str(pd) '% ' Nnull_str],['d' var_str ' > 0% ' Nperm_str]);
        else
            legend(['All data ' Nall_str],['Sig ' var_str ' > 0% ' Nperm_str]);
        end
        ylim([-5*pd 5*pd])
    end
    
    % Plot number of cells
    
    if plot_Sch_sens_cumulative
        Nsigcells = zeros(size(pd_range));
        for i = 1:length(pd_range)
            [q4_boot, pvals1, p0] = calc_significant_permutation([1],pvals_Cave_range(:,:,i),alpha,bad_any,1,do_bh);
            Nsigcells(i) = sum(q4_boot);
        end
        figure; plot(pd_range,Nsigcells / length(q4_boot)*100);
        ylabel(['% cells with ' var_str ' < d' var_str]); xlabel(['d' var_str ' (%)']);
    end
    
    clear s_perm s_boot
    out = vars_push;
    
end