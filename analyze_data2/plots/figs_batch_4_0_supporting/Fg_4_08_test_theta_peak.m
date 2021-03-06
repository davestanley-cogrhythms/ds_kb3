

function [wrkspc_buffer, group] = Fg_4_08_test_theta_peak(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,opts_exclude,perm2pls_split_plusminus);

    % Pls switches
    freqband_stats = [16 20];
    freqband_stats_perm = [16 20];
%     freqband_stats = [1 10];
%     freqband_stats_perm = [1 10];
%     freqband_stats = [10 12];
%     freqband_stats_perm = [9 11];
%     freqband_stats_perm = [8 10];
%     freqband_stats_perm = [10 12];
%     freqband_stats = [30 40];
%     freqband_stats_perm = [30 40];
%     freqband_stats_perm = [80 100];


    % More pls switches
    opts_pls = Opts_Pls;
    opts_pls.plotmode = 1;                   % 1-SFC; 2-PSD; 3-time series; 4-phase; 5-SpkPSD; 6- Abs(Cave*exp(i*phi)); 7 - angle(Cave*exp(i*phi)); 8-Cave*exp(i*phi))
    opts_pls.permdat2pls = 0; % Instead of using Cave from sfc_mode, use Cave1 and Cave2 data from perm_mode
    opts_pls.perm2pls = 1;                   % Instead of using Cave from sfc_mode, show differences between Cave1 and Cave2 normalized by shuffle (zscore)
        opts_pls.perm2pls_do_bh = 0;         % Do bh stepup instead of basic test
        opts_pls.perm2pls_dophi = 0;
        opts_pls.perm2pls_return_mode = 1;                 % Return mode of perm2pls (4=zscore)
        opts_pls.perm2pls_allow_signed = 0;                % If 1, alows perm2pls return to be +ve or -ve; else always positive.
        opts_pls.perm2pls_split_plusminus = perm2pls_split_plusminus;             % 0-Return everything; 1-return pos+negative in separate columns; 2-return only positive (pls.*Ramp(Diff)); 3-return only negative (pls.*Ramp(-Diff)); 4-return only significant cells (doesn't really have a use now)
        opts_pls.sort_pls = 0;               % Sort into preferred and non-preferred
        opts_pls.swap_pls = [];
        %opts_pls.swap_pls = [2,4;5,7];      
        opts_pls.do_diff = 0;                % Take difference between adjacent pls (pls(:,:,2) - pls(:,:,1))

    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 0;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 0;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dcs.stats_mode = 0;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    
    
plot_on_spect = 0;
plot_on_spectrogram = 0;
plot_on_scatter = 0;
plot_on_bargraph = 0;

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


groupmode = 0;
% Use default grouping (all pairs, enumerate over ctgs)
switch groupmode
    case {0}
        group = group0;
        sp = ~bad_any(:);
        
end

% Load data into groups
%group = group.query_legend(group0);
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);

if group_do_merge
    % Flip if ctgsetli mode = 4 or 5.
    temp=1:length(group);
    N=length(group);
    temp=flipud(reshape(1:N,2,floor(N/2)));
    temp=temp(:);
    group = grouppairs_merge(group(temp),opts_pls.perm2pls_dophi,grouppmerge_do_percent);
end

plot_inds = 1:length(group);

% Test Plot groups, at last
opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1};
if plot_on_spect
    
    inds = 1:length(group);
%     inds = [1:4];
%     inds = [5:8];
    figure;
    [h1] = plot_matrix3D_custstruct([],group(plot_inds),opts_PM3D,opts_PM3Dcs);
    
end


% Test spectrogram
if plot_on_spectrogram
    plot_spectrogram_custstruct(group(plot_inds),opts_PM3Dsp);
    
    
end


% Test bargraph

if plot_on_bargraph
    fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group(plot_inds),opts_PSC);
end


% Test Plot scattergroups

if plot_on_scatter
    figure; % subplot(221); 
    ind1 = 1;
    ind2 = 3;
    x = pls_stats(~bad_any,group(ind1).ctgs);
    y = pls_stats(~bad_any,group(ind2).ctgs);
    plott_fit(x,y,'k.');
    hold on; plot_scattergroups(x,y,[group(ind1),group(ind2)],~bad_any);
    xlabel(['FFC ' group(ind1).legend]); ylabel(['FFC ' group(ind2).legend]); % Might need fixing
    
    
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
