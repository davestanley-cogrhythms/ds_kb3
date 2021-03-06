function [wrkspc_buffer ] = Fig_9e_time_series_plots_ctgsetli_mode_7(wrkspc_buffer,sfc_mode,sfc_mode2, ...
    curr_stage_sfc,freqband_stats,opts_exclude,opts_pls,opts_pls2,perm_mode_sp,curr_stage_sp,freqband_stats_perm, ...
    opts_perm,running_mode)



    opts_PM3Dsp.paperfig_mode=1;
    
    

    cell_group_mode = 2;    % 1-enumerate ctgs (obsolete); 2-Grouped deciders; 3-All cells (same as 1 for Fig9b)
    do_AB_symmetry = 0;
    sort_on = 0; % Setting to 0 completely disables all sorting
    do_sort_irr = 0; % 1-process irrelevants separately; 0-sort everything based on relevant category(maybe biased)
    
    do_unitsCI_vs_electBNDRY = 0;
    [wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
    pls_CI = out_pls.pls;
    pls_stats = out_pls.pls_stats;
    abscissa1 = out_pls.abscissa;
    bad_any1 = out_pls.bad_any;
    funames1D1 = out_pls.funames1D;
    mypairs1 = out_pls.mypairs;
    group1 = out_pls.group;
    
    % % Load PLS for boundary response spectrogram
    [wrkspc_buffer, out_pls2] = load_pls(wrkspc_buffer,sfc_mode2,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls2);
    if opts_pls2.permdat2pls factor = -1; else factor = 1; end
    pls2 = factor * out_pls2.pls;
    pls_stats2 = factor* out_pls2.pls_stats;
    abscissa21 = out_pls2.abscissa;
    abscissa22 = out_pls2.abscissa2;
    
    group2 = out_pls2.group;
    bad_any2 = out_pls2.bad_any;
    
  
    %%
    % Group1 stuff
    grt1 = group1(1);
    grt1.criteria = [2]; grt1.criteria_alt = []; grt1.criteria_sfc = [];
    grt1.xdata = group1(1).xdata;
    grt1.data_name = ['Category Index'];
    grt1.xlims_desired = [-1.4 1.5];
    
    sp = true(size(pls_CI,2),1);
    
    grt1 = get_grouped_cells_and_data(grt1,sp,pls_CI,abscissa1,mypairs1,bad_any1,opts_pls.plotmode,freqband_stats,funames1D1);
    
    %
    % Group2 stuff
    grt2 = group2(1);
    grt2.criteria = [2]; grt2.criteria_alt = []; grt2.criteria_sfc = [];
    grt2.xdata = group1(1).xdata;
    grt2.data_name = ['Boundary response'];
    grt2.xlims_desired = [-1.4 1.5];
    grt2.ylims_desired = [0 50];
    %grt2.zlims_desired = [-0.25 0.25];
    
    sp2 = true(size(pls2,2),1);
    
    grt2 = get_grouped_cells_and_data(grt2,sp2,pls2,out_pls2.abscissa,out_pls2.mypairs,out_pls2.bad_any,opts_pls.plotmode,freqband_stats,out_pls2.funames1D,out_pls2.abscissa2);
    grt2.data = pls2(:,7,1,:);
    
    
    
    
    
    
    % Map FFC to units
    ulpairs = build_unit_LFP_pairs_all(wrkspc_buffer.currmd.md);
    unit2lfp = ulpairs(2,:);
    mypairs_units = out_pls.mypairs(:,2);
    mypairs_FFC = out_pls2.mypairs;
    
	%% Do stuff
    
    if running_mode == 1            % Plot  
        %%  Plot
        figure;
        set(gcf,'PaperPositionMode','auto');

        fprintf('Adjust figure as needed.');

        pause
        set(gcf,'Visible','off');
        currfname = '';
        currfigname = 'Fig9d';
        mydate = datestr(datenum(date),'yy/mm/dd'); mydate = strrep(mydate,'/','');
        c=clock;
        sp = ['d' mydate '_t' num2str(c(4),'%10.2d') '' num2str(c(5),'%10.2d') '' num2str(round(c(6)),'%10.2d')];
        sp = [sp '__' currfname '_' currfigname];

        basepath = '~/figs_tosave';
        mkdir(fullfile(basepath,sp));

        for i = 1:length(mypairs_FFC)
            %warning('Identify low firing rate cells');
            if ~bad_any2(i)
                u1 = find(unit2lfp == mypairs_FFC(i,1));
                u2 = find(unit2lfp == mypairs_FFC(i,2));

                subplot(211);
                grt2.data = pls2(:,i,1,:);
                plot_spectrogram_custstruct(grt2,opts_PM3Dsp);
                title(['Mypairs ' num2str(mypairs_FFC(i,1)) ',' num2str(mypairs_FFC(i,2))]);

                subplot(212);
                u = [u1, u2];
                clear grt1_all
                if ~isempty(u)
                    for j = 1:length(u)
                        grt1_all(j) = grt1;
                        grt1_all(j).data = pls_CI(:,u(j),1);
                        grt1_all(j).legend = ['Unit ' num2str(u(j))];
                    end

                    [~, out.PM3Dcs] = plot_matrix3D_custstruct([],[grt1_all]);
                    colorbar
                end

                tic
                print(gcf,'-dpng',fullfile(basepath,sp,['fig' num2str(i)]))
                %print(gcf,fullfile(basepath,sp,['fig' num2str(i)]))
                toc
                %pause
                clf;

            end
        end
    elseif running_mode == 2            % Analyze
    
        %% correlation coefficient

        ctg1 = 6;
        ctg2 = 3;
        corrcoefs = [];
        plot_debug = 0;

        % Resizing parameters
        tstart = 0.0;
        tstop = 1.5;
        ind = abscissa22 >= tstart & abscissa22 < tstop;

        for i = 1:length(mypairs_FFC)
            i;
            %warning('Identify low firing rate cells');
            if ~bad_any2(i)
                u1 = find(unit2lfp == mypairs_FFC(i,1));
                u2 = find(unit2lfp == mypairs_FFC(i,2));

                d1 = pls2(:,i,ctg2,:);
                d1 = squeeze(d1);

                % Resize
                d1 = d1(:,ind);


                u = [u1, u2];
                if ~isempty(u)
                    corrcoef_temp2 = [];
                    for j = 1:length(u)
                        j;
                        d2 = pls_CI(:,u(j),ctg1);

                        % Interpolate time series
                        d2ds = interp1(abscissa1(:), d2,abscissa22);

                        % Resize
                        d2ds = d2ds(ind);

                        if plot_debug
                            %%
                            figure; 
                            plot(abscissa1,d2);
                            hold on; plot(abscissa22(ind),d2ds,'r');
                            legend('Original','Downsampled')
                            pause
                        end

                        corrcoef_temp = zeros(size(d1,1),1);
                        for k = 1:size(d1,1)
                            [R] = corrcoef(d1(k,:),d2ds);
                            corrcoef_temp(k) = R(2);
                        end
                        
                        corrcoef_temp2 = [corrcoef_temp2, corrcoef_temp];       % Inefficient...
                        

                    end
                    
                    corrcoefs = [corrcoefs, mean(corrcoef_temp,2)];     % Average all unit responses for each pair
                    %corrcoefs = [corrcoefs, corrcoef_temp];             % Include each unit separately
                    

                end


            end
        end
        
        figure; plott_matrix3D(abscissa21,corrcoefs)
    end
    
    %%
    
%     corrcoefs_100 = corrcoefs;
%     corrcoefs_60 = corrcoefs;
%     figure; 
%       subplot(211); 
%       plott_matrix3D(abscissa21,cat(3,corrcoefs_100,corrcoefs_60));
%       title('SchA 50%'); legend('60% ctgs z-score','100% ctgs z-score');
%      for i = 1:52; p(i) = ranksum(corrcoefs_100(i,:), corrcoefs_60(i,:)); end
%     subplot(212);
%       plot(abscissa21,p);
%     
   
    
    
   
end

