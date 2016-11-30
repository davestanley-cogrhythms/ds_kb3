function [wrkspc_buffer, group] = Fg_4_01_alpha_pls(wrkspc_buffer,sfc_mode, perm_mode,opts_exclude, opts_perm,freqband_stats,freqband_stats_perm,groupmode)

%%

% Setup params 
run_setdefaultfig
addpath(genpath('./funcs_plots_preferred/'));
addpath(genpath('./funcs_generalized/'));
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end

% SFC mode
% sfc_mode = 22.4018111;
% perm_mode = 22.4018111;
% sfc_mode =  2.3515111;
% perm_mode = 2.3015111;
% perm_mode = 5;
% sfc_mode = 41.6514111;
% perm_mode = 41.6514111;
% sfc_mode = 52.7002010;
% perm_mode = 52.7000010;


% Stage selection
curr_stage_sfc = 3;
curr_stage_sp = 3;

    % Pls switches
%     freqband_stats = [16 20];
%     freqband_stats_perm = [16 20];
%     freqband_stats = [1 10];
%     freqband_stats_perm = [1 10];
%     freqband_stats = [10 12];
%     freqband_stats_perm = [9 11];
%     freqband_stats_perm = [8 10];
%     freqband_stats_perm = [10 12];
%     freqband_stats = [30 40];
%     freqband_stats_perm = [30 40];
%     freqband_stats_perm = [80 100];

    % Set up pls fv (popts)
    % Unit exclusion
%     opts_exclude = Opts_Exclude;
%     opts_exclude.exclude_clipping = 1;
%     opts_exclude.exclude_60 = 0;
%     opts_exclude.exclude_nans = 1;
%     opts_exclude.excludeL = 1;
%     opts_exclude.excludeO = 0; 
%     opts_exclude.remove_dependent = 0;       % Remove dependent electrode pairs for FFC


    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 1; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 0;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 4;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))

    % Permutation test options
%     opts_perm = Opts_Perm;
%     opts_perm.do_bh0 = 0;
%     opts_perm.do_phi = 0;
%     opts_perm.split_plusminus = 1;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dcs.stats_mode = 0;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    
    warning('To do - compare sch sensitive electrodes ability to attenuate irrelevant ctg response to other electrodes.');

% Plot switches
plot_on_spect = 1;
plot_on_spectrogram = 0;
plot_on_scatter = 0;
plot_on_bargraph = 0;
plot_on_histogram = 0;
plot_on_real_vs_imag = 0;


group_do_merge = 0;
    grouppmerge_do_percent = 1;


        
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

% pls = pls * -1;
% pls_stats = pls_stats * -1;


% Build group

% Use default grouping (all pairs, enumerate over ctgs)
group = group0;


% Load sp's
[wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm,opts_exclude);
sp = out_perm.sig_cells;


% Load bads perm
[bad_any_perm] = load_bads(perm_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,out_perm.mypairs,out_perm.funames1D);

% Map sp's as needed
[sp] = map_sp(perm_mode, sfc_mode,out_perm.mypairs,mypairs,sp,wrkspc_buffer.currmd.md,bad_any_perm,bad_any);

% Create group template
mycrit = [2*ones(1,size(sp,2))];
grt = group0(1);
grt.criteria = mycrit; grt.criteria_alt = []; grt.criteria_sfc = []; grt.ctgs = 1;

% Run a simple test
clear group
switch groupmode
    case 0      % Plot all schemes
        i=0;
        % Sch A&B
        i=i+1; group(i)=grt; group(i).criteria(1)=[1]; group(i).ctgs=1;   % Ctg1-2 deciders
        i=i+1; group(i)=grt; group(i).criteria(2)=[1]; group(i).ctgs=2;   % Ctg1-2 non-deciders
        i=i+1; group(i)=grt; group(i).criteria(3)=[1]; group(i).ctgs=3;   % Ctg1-2 deciders
        i=i+1; group(i)=grt; group(i).criteria(4)=[1]; group(i).ctgs=4;   % Ctg1-2 non-deciders
        group = group.query_legend(group0);
        
    case 1      % Plot just scheme A (preferred & non-preferred)
            i=0;
            % Sch A
            i=i+1; group(i)=grt; group(i).criteria(1:2)=[1,0]; group(i).ctgs=9;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:2)=[0,0]; group(i).ctgs=[9];   % Ctg1-2 non-deciders
            i=i+1; group(i)=grt; group(i).criteria(1:2)=[0,1]; group(i).ctgs=9;   % Ctg1-2 deciders
            group = group.query_legend(group0);
            
            i=0;
            i=i+1; group(i).legend = 'Pref Cat';
            i=i+1; group(i).legend = 'No Pref';
            i=i+1; group(i).legend = 'Pref Dog';
    
        
    case 2
            i=0;
            % Sch A
            i=i+1; group(i)=grt; group(i).criteria(3:4)=[1,0]; group(i).ctgs=10;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(3:4)=[0,0]; group(i).ctgs=[10];   % Ctg1-2 non-deciders
            i=i+1; group(i)=grt; group(i).criteria(3:4)=[0,1]; group(i).ctgs=10;   % Ctg1-2 deciders
            group = group.query_legend(group0);
            
            i=0;
            i=i+1; group(i).legend = 'Pref Fat';
            i=i+1; group(i).legend = 'No Pref';
            i=i+1; group(i).legend = 'Pref Thin';
              
    case 3
            i=0;
            % Sch A Cat
            i=i+1; group(i)=grt; group(i).criteria(1:2)=[1,0]; group(i).ctgs=1;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:2)=[1,0]; group(i).ctgs=2;   % Ctg1-2 non-deciders
%             i=i+1; group(i)=grt; group(i).criteria(1:2)=[0,1]; group(i).ctgs=2;   % Ctg1-2 deciders
%             i=i+1; group(i)=grt; group(i).criteria(1:2)=[0,1]; group(i).ctgs=1;   % Ctg1-2 non-deciders
            
            i=0;
            i=i+1; group(i).legend = 'Pref Cat:Cat';
            i=i+1; group(i).legend = 'Pref Cat:Dog';
%             i=i+1; group(i).legend = 'Pref Dog:Dog';
%             i=i+1; group(i).legend = 'Pref Dog:Cat';


    case 4
            i=0;
            i=i+1; group(i)=grt; group(i).criteria(1:2)=[0,1]; group(i).ctgs=2;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(1:2)=[0,1]; group(i).ctgs=1;   % Ctg1-2 non-deciders
            group = group.query_legend(group0);
            i=0;
            i=i+1; group(i).legend = 'Pref Dog:Dog';
            i=i+1; group(i).legend = 'Pref Dog:Cat';
            
    case 5
            i=0;
            % Sch B
            i=i+1; group(i)=grt; group(i).criteria(3:4)=[1,0]; group(i).ctgs=3;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(3:4)=[1,0]; group(i).ctgs=4;   % Ctg1-2 non-deciders
%             i=i+1; group(i)=grt; group(i).criteria(3:4)=[0,1]; group(i).ctgs=4;   % Ctg1-2 deciders
%             i=i+1; group(i)=grt; group(i).criteria(3:4)=[0,1]; group(i).ctgs=3;   % Ctg1-2 non-deciders
            i=0;
            i=i+1; group(i).legend = 'Pref Fat:Fat';
            i=i+1; group(i).legend = 'Pref Fat:Thin';
%             i=i+1; group(i).legend = 'Pref Thin:Thin';
%             i=i+1; group(i).legend = 'Pref Thin:Fat';

                        
    case 6
            i=0;    
            i=i+1; group(i)=grt; group(i).criteria(3:4)=[0,1]; group(i).ctgs=4;   % Ctg1-2 deciders
            i=i+1; group(i)=grt; group(i).criteria(3:4)=[0,1]; group(i).ctgs=3;   % Ctg1-2 non-deciders
            group = group.query_legend(group0);
            i=0;
            i=i+1; group(i).legend = 'Pref Thin:Thin';
            i=i+1; group(i).legend = 'Pref Thin:Fat';
            
end



% Load data into groups
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);

if group_do_merge
    % Flip if ctgsetli mode = 4 or 5.
    temp=1:length(group);
    N=length(group);
    temp=flipud(reshape(1:N,2,floor(N/2)));
    temp=temp(:);
    group = grouppairs_merge(group(temp),opts_pls.perm2pls_dophi,grouppmerge_do_percent);
end

% Test Plot groups, at last
opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
if plot_on_spect
    
    switch groupmode
        case {1,2,3,4,5,6}
            inds = 1:length(group);
            figure; [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
        case {}
            Nhalf=floor(length(group)/2);
            figure;[h1] = plot_matrix3D_custstruct([],group(1:Nhalf),opts_PM3D,opts_PM3Dcs);
            figure; [h1] = plot_matrix3D_custstruct([],group(Nhalf+1:end),opts_PM3D,opts_PM3Dcs);
    end
        
    
end


% Test spectrogram
if plot_on_spectrogram
    plot_spectrogram_custstruct(group(:),opts_PM3Dsp);
end


% Test bargraph

if plot_on_bargraph
%     hmask = true(length(group));
%     hmask(logical(eye(size(hmask)))) = false;
%     opts_PSC.hmask = hmask;
%     opts_PSC.multicomparisons_mode = 2;     % 2-default (holm-bonferroni)
%     opts_PSC.multicomparisons_mode = 0; warning('Multiple comparisons turned off');
    fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
end


% Test Plot scattergroups

if plot_on_scatter
    figure; % subplot(221); 
    ind1 = 1;
    ind2 = 2;
    x = pls_stats(~bad_any,group(ind1).ctgs);
    y = pls_stats(~bad_any,group(ind2).ctgs);
    plott_fit(x,y,'k.');
    hold on; plot_scattergroups(x,y,[group(ind1),group(ind2)],~bad_any);
    xlabel(['FFC ' group(ind1).legend]); ylabel(['FFC ' group(ind2).legend]); % Might need fixing
    hold on; plot([0,1],[0,1],'k','LineWidth',5);
    
    fprintf(['L1 Slope SchA Rel Vs Irrel= ' num2str(std(pls_stats(~bad_any,3))/std(pls_stats(~bad_any,1))) '\n']);
    fprintf(['L1 Slope SchB Rel Vs Irrel= ' num2str(std(pls_stats(~bad_any,4))/std(pls_stats(~bad_any,2))) '\n']);
    
    
end

do_pca = 0;
if do_pca
   [coeff,score,latent] = princomp([x(:),y(:)]);
    ax1 = coeff(:,1);% *sqrt(latent(1));
    ax2 = coeff(:,2);% *sqrt(latent(2));
    hold on; plot([0 coeff(1,1)],[0 coeff(2,1)],'r','LineWidth',2);
    hold on; plot([0 coeff(1,2)],[0 coeff(2,2)],'r','LineWidth',2); 
end

% Plot histogram
if plot_on_histogram
    %plot_histplot(group,'plot_type','hist_paired');
    plot_histplot(group,'plot_type','hist_autogrouped');
end

% Plot real vs imag

if plot_on_real_vs_imag
    valid_plot_types = {'real_imag_scatter','angle_cdf','angle_pdf','mag_vs_angle_scatter'};
    h = plot_phi_analysis(group,'plot_type',valid_plot_types{4});
end


end