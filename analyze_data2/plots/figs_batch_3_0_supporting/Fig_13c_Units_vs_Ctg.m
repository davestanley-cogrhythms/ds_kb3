function [wrkspc_buffer, group] = Fig_13c_Units_vs_Ctg(wrkspc_buffer,do_units,do_percent_change,do_ABsymmetry,curr_stage_sfc,curr_stage_sp)

% 
% do_units=1;             % 0-FFC; 1-Units time series; 2-units by stages
% do_percent_change = 1;
% do_ABsymmetry = 1;
% 
% 
% % Stage selection
% curr_stage_sfc = 3;
% curr_stage_sp = 3;

% These are 2 different modes; Enable one and disable the other
does_ensemble_confer_stronger_scheme_preference = 1;    
does_ensemble_confer_stronger_unit_response = 0;


% Plot switches
plot_on_spect = 1;
plot_on_bargraph = 1;

if do_units == 0            % FFC
    % SFC mode FFC
    pername = 'Units';
    sfc_mode = 22.401411101;
    perm_mode = 52.700001001;    

elseif do_units == 1        % Units time series
    % SFC mode units
    pername = 'Ens.';
    sfc_mode = 52.700201001;
    perm_mode = 22.401411101;
    plot_on_bargraph = 0;
    
elseif do_units == 2        % Units by stages
    pername = 'Ens.';
    sfc_mode = 52.700001001;
    perm_mode = 22.401411101;
    plot_on_spect = 0;
end

if do_units == 1
    curr_stage_sfc=4;
end

    % Pls switches
    freqband_stats = [16 20];
    freqband_stats_perm = [16 20];
%     freqband_stats = [1 10];
%     freqband_stats_perm = [1 10];
%     freqband_stats = [10 12];
%     freqband_stats_perm = [10 12];
%     freqband_stats = [30 40];
%     freqband_stats_perm = [30 40];


    % Set up pls fv (popts)
    % Unit exclusion
    opts_exclude = Opts_Exclude;
    opts_exclude.exclude_clipping = 1;
    opts_exclude.exclude_60 = 0;
    opts_exclude.exclude_nans = 1;
    opts_exclude.excludeL = 0; 
    opts_exclude.excludeO = 0; 
    opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    
    
    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))
        
    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 1;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    opts_perm.alpha_bh0 = 0.2;
%     opts_perm.alpha_bh0 = 0.05;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dcs.stats_mode = 0;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;


    % Load pls
    [wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode,curr_stage_sfc,freqband_stats,opts_exclude,opts_pls);
    pls = out_pls.pls;
    pls_stats = out_pls.pls_stats;
    abscissa = out_pls.abscissa;
    abscissa2 = out_pls.abscissa2;
    bad_any = out_pls.bad_any;
    funames1D = out_pls.funames1D;
    mypairs = out_pls.mypairs;
    group0 = out_pls.group;
    clear out_pls

    % Load sp's
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm,opts_exclude);
    sp = out_perm.sig_cells;

    % Load bads perm
    [bad_any_perm] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);
    %
    % Map sp's as needed
    [sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any);

    % Create group template
    mycrit = [2*ones(1,size(sp,2))];
    grt = group0(1);
    grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;

    % Run a simple test
    clear group
    i=0;
    if does_ensemble_confer_stronger_scheme_preference
        i=i+1; group(i)=grt; group(i).criteria(1:2)=[1 0]; group(i).ctgs=1; group(i).legend = ['$\rm \bf C/D \thinspace ' pername '  Pref $']; 
        i=i+1; group(i)=grt; group(i).criteria(1:2)=[1 0]; group(i).ctgs=2; group(i).legend = ['$\rm \bf C/D \thinspace ' pername '  NPref$'];
    end
    if does_ensemble_confer_stronger_unit_response
        i=i+1; group(i)=grt; group(i).criteria(1:2)=[1 2]; group(i).ctgs=1; group(i).legend = ['$\rm \bf C/D \thinspace ' pername '  $']; 
        i=i+1; group(i)=grt; group(i).criteria(1:2)=[0 2]; group(i).ctgs=1; group(i).legend = ['$\rm \bf \overline{C/D} \thinspace ' pername ' $'];
    end
    
    if does_ensemble_confer_stronger_scheme_preference
        i=i+1; group(i)=grt; group(i).criteria(1:2)=[0 1]; group(i).ctgs=2; group(i).legend = ['$\rm \bf F/T \thinspace ' pername '  Pref $']; 
        i=i+1; group(i)=grt; group(i).criteria(1:2)=[0 1]; group(i).ctgs=1; group(i).legend = ['$\rm \bf F/T \thinspace ' pername '  Npref $']; 
    end
    if does_ensemble_confer_stronger_unit_response
        i=i+1; group(i)=grt; group(i).criteria(1:2)=[2 1]; group(i).ctgs=2; group(i).legend = ['$\rm \bf F/T \thinspace ' pername '  $']; 
        i=i+1; group(i)=grt; group(i).criteria(1:2)=[2 0]; group(i).ctgs=2; group(i).legend = ['$\rm \bf \overline{F/T} \thinspace ' pername ' $'];
    end
    
    for i = 1:length(group); group(i).interpreter = 'latex'; end
    
    if do_ABsymmetry
        clear group2
        N=floor(length(group)/2);
        for i = 1:N
            group2(i) = group_merge(group(i),group(i+N));
        end
        
        group_unmerged=group;    % Save original group
        group=group2;
        
        i=0;
        if does_ensemble_confer_stronger_scheme_preference
            i=i+1; group(i).legend = 'Ens. Pref';
            i=i+1; group(i).legend = 'Ens. Non-pref';
        end
        if does_ensemble_confer_stronger_unit_response
            i=i+1; group(i).legend = 'Ens.';
            i=i+1; group(i).legend = 'Non-Ens.';
        end
        
    end
    
    % Load data into groups
    group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);
    
    % Test Plot groups, at last
    opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
    if plot_on_spect
        N = length(group);
        
        inds = 1:N;
        figure; [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
        
%         inds = 1:round(N/2);
%         figure; [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
%         
%         inds = round(N/2)+1:N;
%         figure; [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    end

    % Test bargraph
    if plot_on_bargraph
        fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
    end


end