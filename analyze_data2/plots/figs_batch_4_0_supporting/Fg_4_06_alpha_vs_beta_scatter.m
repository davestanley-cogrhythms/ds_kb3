
%% Run!
% Use this code for copying and pasting into supporting functions


function [wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter(wrkspc_buffer,sfc_mode1,sfc_mode2,curr_stage1,curr_stage2,freqband_stats1,freqband_stats2,plotmode1,plotmode2,opts_pls,opts_exclude)



    %
    % FFC
    sfc_mode = sfc_mode1;
    perm_mode = sfc_mode1;
    freqband_stats = freqband_stats1;
    freqband_stats_perm = freqband_stats1;
    curr_stage_sfc = curr_stage1;
    curr_stage_sp = curr_stage1;
    plotmode = plotmode1;
    
    opts_pls.target_pls_format = 0; % Convert pls to match this format!

    [wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter_supporting(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,plotmode,opts_pls,opts_exclude);
    pls_stats_ffc = pls_stats;
    group_ffc = group;
    bad_ffc = bad_any;

    %
    % PSD
    sfc_mode = sfc_mode2;
    perm_mode = sfc_mode2;
    freqband_stats = freqband_stats2;
    freqband_stats_perm = freqband_stats2;
    curr_stage_sfc = curr_stage2;
    curr_stage_sp = curr_stage2;
    plotmode = plotmode2;
    
    if sfc_mode1 ~= sfc_mode2
        opts_pls.target_pls_format = sfc_mode1; % Convert pls to match this format!
    else
        opts_pls.target_pls_format = 0;
    end

    [wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter_supporting(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,plotmode,opts_pls,opts_exclude);
    pls_stats_psd = pls_stats;
    group_psd = group;
    bad_psd = bad_any;


    if any(bad_ffc ~= bad_psd); error('Bads mismatch'); end
    figure; % subplot(221);
    indx = 9;
    indy = 19;
    pls_stats = cat(2,pls_stats_ffc,pls_stats_psd);
    plot_scatterpls_and_groups(pls_stats,[group_ffc],bad_ffc,indx,indy)
    %xlabel('PSD');ylabel('FFC');

end


function [wrkspc_buffer, pls_stats, group, group0, bad_any] = Fg_4_06_alpha_vs_beta_scatter_supporting(wrkspc_buffer,sfc_mode,perm_mode,curr_stage_sfc,curr_stage_sp,freqband_stats,freqband_stats_perm,plotmode,opts_pls,opts_exclude)


    

%     freqband_stats = [10 12];
%     freqband_stats_perm = [9 11];
%     freqband_stats_perm = [8 10];
%     freqband_stats_perm = [10 12];
%     freqband_stats = [30 40];
%     freqband_stats_perm = [30 40];
%     freqband_stats_perm = [80 100];


    opts_pls.plotmode = plotmode;


    % Permutation test options
    opts_perm = Opts_Perm;
    opts_perm.do_bh0 = 0;
    opts_perm.do_phi = 0;
    opts_perm.split_plusminus = 1;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dcs.stats_mode = 0;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    
    
% Plot switches
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
sfc_mode = out_pls.sfc_mode;
clear out_pls

% pls = pls * -1;
% pls_stats = pls_stats * -1;


% Build group
groupmode = 1;

% Use default grouping (all pairs, enumerate over ctgs)
switch groupmode
    case 0
        group = group0;
        sp = ~bad_any(:);
    case 1
    % Load sp's
    [wrkspc_buffer, out_perm] = load_pr(wrkspc_buffer,perm_mode,curr_stage_sp,freqband_stats_perm,bad_any,opts_perm);
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
    i=0;
    i=i+1; group(i)=grt; group(i).criteria(1:4)=[1 0 0 0]; group(i).ctgs=1;   % 
    i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 1 0 0]; group(i).ctgs=1;   % 
    i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 0 1 0]; group(i).ctgs=1;   % 
    i=i+1; group(i)=grt; group(i).criteria(1:4)=[0 0 0 1]; group(i).ctgs=1;   % 
    

end

% Load data into groups
group = group.query_legend(group0);
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
    
    
    inds = 1:length(group);
    inds = [1:4];
%     inds = [5:8];
    figure;
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    
end


% Test spectrogram
if plot_on_spectrogram
    plot_spectrogram_custstruct(group(:),opts_PM3Dsp);
    
    
end


% Test bargraph

if plot_on_bargraph
    fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
end


% Test Plot scattergroups
if plot_on_scatter
    figure; % subplot(221); 
    ind1 = 1;
    ind2 = 3;
    plot_scatterpls_and_groups(pls_stats,group,bad_any,ind1,ind2)
end

do_pca = 0;
if do_pca
   [coeff,score,latent] = princomp([x(:),y(:)]);
    ax1 = coeff(:,1);% *sqrt(latent(1));
    ax2 = coeff(:,2);% *sqrt(latent(2));
    hold on; plot([0 coeff(1,1)],[0 coeff(2,1)],'r','LineWidth',2);
    hold on; plot([0 coeff(1,2)],[0 coeff(2,2)],'r','LineWidth',2); 
end

end
