
%% All figures

run_setdefaultfig
addpath(genpath('./funcs_plots_preferred/'));
fv.mn = {'L','O'};
fv.exclude_clipping = 1;
fv.exclude_60 = 0;
fv.exclude_lowfiring = 1;
fv.exclude_nans = 1;

fv.plot_SFC_spectra = 0;
fv.compare_boundary_AvB = 0;
fv.plot_SFC_spectra_stats = 0;
fv.plot_Sch_sens = 0;
fv.plot_spectrogram = 0;        % Make sure this is zero, since they can mess with our badfiles settings if set otherwise
fv.plot_spectrogram_vs_ev=0;    % This too.

% Figure plotting structs
paperfig_mode=1;
fv.opts_PM3Dcs.paperfig_mode=paperfig_mode;
fv.opts_PM3Dcs.do_sgolay=0;
fv.opts_PM3Dcs.remove_dependent=0;
fv.opts_PSC.paperfig_mode=paperfig_mode;
fv.opts_PSC.remove_dependent=0;
fv.opts_PPS.paperfig_mode=paperfig_mode;
fv.opts_PPS.hide_boot_results=1;
fv.opts_PPS.sens_chosen_condition=1;

%fv.freqband_stats_hbeta = [19 39];
fv.freqband_stats_hbeta = [20 32];
fv.freqband_stats_alpha = [10 12];
fv.freqband_stats_lbeta = [16 20];


fv.swap_in_PEV = 1;        % Use percent explained variance instead of sensitivity measures
fv.preferred_mode = 3;     % 1 for statistics (ttest, not working in plots_preferred_unitpairs.m); 2 for raw threshold; 3 for quantile


fv.plotmode = 1;
fv.permdat2pls = 0;            % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
fv.perm2pls = 0;               % Instead of showing raw SFC values, show % cells successfully passing permutation test
        fv.perm2pls_do_bh = 0;        % Do bh stepup instead of basic test
        fv.perm2pls_dophi = 0;
        fv.sort_pls = 0;           % Sort into preferred and non-preferred
        fv.do_diff = 0;            % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
        
fv.savetype='png';
fv.savenames= {'fig1','fig2','fig3','fig4','fig5','fig6','fig7','fig8','fig9','fig10','fig11','fig12','fig13','fig14'};
% Turn off all script plotting