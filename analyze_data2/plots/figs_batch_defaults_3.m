
%% All figures

function fv = figs_batch_defaults_3

    run_setdefaultfig

    % General input vars defaults
    fv.sfc_mode = 22.4014111;
    fv.perm_mode = fv.sfc_mode;
    
    fv.curr_stage_sfc = 2;
    fv.curr_stage_sp = 2;
    
    fv.freqband_stats = [16 20];

    % Pack option structures
    fv.opts_exclude = Opts_Exclude;
    fv.opts_pls = Opts_Pls;
    fv.opts_perm = Opts_Perm;

    
    % Figure plotting structs
    paperfig_mode=1;
    fv.opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
    fv.opts_PM3Dcs.paperfig_mode=paperfig_mode;
    fv.opts_PM3Dcs.do_sgolay=0;
    fv.opts_PM3Dcs.remove_dependent=0;
    fv.stats_mode.remove_dependent=0;
    fv.opts_PSC.paperfig_mode = paperfig_mode;
    fv.opts_PSC.remove_dependent = 0;
    fv.opts_PSC.do_stats = 1;
    fv.opts_PSC.do_diagtest = 0;
    fv.opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    fv.opts_PSC.opts_diagtest.threshold = 1.0; 
    
    
    
    
    fv.opts_PPS.paperfig_mode=paperfig_mode;
    fv.opts_PPS.hide_boot_results=1;
    fv.opts_PPS.sens_chosen_condition=1;

    % Plot switches
    fv.plot_SFC_spectra = 0;
    fv.plot_Sch_sens = 0;
    
    %Group switches
    fv.group_by_perm = 1;   %0-Use default grouping (all pairs, enumerate over ctgs); 1-Group based on permutation test output.


    % Figure save settings        
    fv.savetype='png';
    fv.savenames= {'fig1','fig2','fig3','fig4','fig5','fig6','fig7','fig8','fig9','fig10','fig11','fig12','fig13','fig14'};
    % Turn off all script plotting

end