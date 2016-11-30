function [wrkspc_buffer, group] = Fig_17_FFC_partial_coherence_index(wrkspc_buffer, sfc_mode1, sfc_mode2,curr_stage_sfc1, curr_stage_sfc2,opts_exclude,opts_pls)



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


    % Map sp parameters
    sp_threshold = 10;
    
    % Group options
    do_group_collapse_pls2days = 0;
    
    % Plotting options
    paperfig_mode = 1;
    opts_PM3Dcs.paperfig_mode=paperfig_mode;
    opts_PM3Dcs.stats_mode = 0;
    opts_PM3Dsp.paperfig_mode=paperfig_mode;
    opts_PSC.paperfig_mode = paperfig_mode;
    opts_PSC.remove_dependent = 0;
    opts_PSC.do_stats = 0;
    opts_PSC.do_diagtest = 1;
    opts_PSC.opts_diagtest.testmode = 2;        % 1-Binomial test against expected N; 2-signrank test against thresh
    opts_PSC.opts_diagtest.threshold = 1.0; 
    
%% Plot switches
plot_on_spect = 0;
plot_on_spectrogram = 0;
plot_on_scatter = 0;
plot_on_bargraph = 0;

groupmode = 0;   % 0-Use default grouping (all pairs, enumerate over ctgs);
                 % 1-Group based on permutation test output for diff.
                 % 2-Group based on permutation test output all.
    swap_mode = 0;

group_do_merge = 0;
    grouppmerge_do_percent = 1;


%% Import data if running as a function; define dependent parameters

running_as_script = ~is_calledby;
if ~running_as_script
    pop_workspace(false);
end

        
%% Load pls1
[wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode1,curr_stage_sfc1,freqband_stats,opts_exclude,opts_pls);
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

% Invert pls if needed
[~, sfc_subgroups] = decode_sfc_mode(sfc_mode);
[~, ~, ctgsetli_mode] = build_sfcmode(sfc_mode, sfc_subgroups);
if any(ctgsetli_mode == [4,5]) && opts_pls.do_diff
    pls = pls * -1;
    pls_stats = pls_stats * -1;
end

pls1 = pls;
pls_stats1 = pls;

%% Load pls2
[wrkspc_buffer, out_pls] = load_pls(wrkspc_buffer,sfc_mode2,curr_stage_sfc2,freqband_stats,opts_exclude,opts_pls);
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

% Invert pls if needed
[~, sfc_subgroups] = decode_sfc_mode(sfc_mode);
[~, ~, ctgsetli_mode] = build_sfcmode(sfc_mode, sfc_subgroups);
if any(ctgsetli_mode == [4,5]) && opts_pls.do_diff
    pls = pls * -1;
    pls_stats = pls_stats * -1;
end

pls2 = pls;
pls_stats2 = pls;

%% Merge pls groups
% for pls
pls_cat_select = abs(cat(3,pls1(:,:,1) - pls1(:,:,2),pls1(:,:,3) - pls1(:,:,4)));
pls_bnd_select = -1*(cat(3,pls2(:,:,1) - pls2(:,:,2),pls2(:,:,3) - pls2(:,:,4)));    % -1* 100% - 50%
pls = (pls_bnd_select - pls_cat_select) ./ (abs(pls_bnd_select) + pls_cat_select);

% for pls_stats
pls_stats_cat_select = abs(cat(3,pls_stats1(:,1) - pls_stats1(:,2),pls_stats1(:,3) - pls_stats1(:,4)));
pls_stats_bnd_select = -1*(cat(3,pls_stats2(:,1) - pls_stats2(:,2),pls_stats2(:,3) - pls_stats2(:,4)));
pls_stats = (pls_stats_bnd_select - pls_stats_cat_select) ./ (abs(pls_stats_bnd_select) + pls_stats_cat_select);

%% Build group

% Group based on perm_mode
if ~exist('group','var')
    switch groupmode
        case 0
            % Use default grouping (all pairs, enumerate over ctgs)
            group = group0(1:2);
            sp = ~bad_any(:);

    end
end

%% Load data into groups
group = get_grouped_cells_and_data(group,sp,pls,abscissa,mypairs,bad_any,opts_pls.plotmode,freqband_stats,funames1D,abscissa2);

if do_group_collapse_pls2days
    group = group_collapse_pls2days(group);
end

if group_do_merge
    % Flip if ctgsetli mode = 4 or 5.
    temp=1:length(group);
    N=length(group);
    temp=flipud(reshape(1:N,2,floor(N/2)));
    temp=temp(:);
    group = grouppairs_merge(group(temp),opts_pls.perm2pls_dophi,grouppmerge_do_percent);
end

%% Test Plot groups, at last
opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',0};
if plot_on_spect
    
    do_split = 0;

    N=length(group);
    if get_iscromer && groupmode == 0
        if length(group) > 5
            inds=[1:4,9,10];
        else
            inds=[1,2,5];
        end
    else
        inds = 1:N;
        %inds = 1:2;
%         inds = [1:2,5];
    end
            
    
    
    if N <= 5 || ~do_split
        figure;
        [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    else
        %%
        inds = 1:floor(N/2); figure; [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
        inds = floor(N/2)+1:N; figure; [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D,opts_PM3Dcs);
    end
    
end

if 0
    %%
    figure;
    inds = [1,2,3,4,9,10];
    [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D);
%     figure;
%     inds = 3:4;
%     [h1] = plot_matrix3D_custstruct([],group(inds),opts_PM3D);

end

%% Test spectrogram
if plot_on_spectrogram
    plot_spectrogram_custstruct(group(:),opts_PM3Dsp);
    
    
end


%% Test bargraph

if plot_on_bargraph
    fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);
end


%% Test Plot scattergroups

if plot_on_scatter

    figure; % subplot(221); 
    indx = 1;
    indy = 5;
    plot_scatterpls_and_groups(pls_stats,group([1,5]),bad_any,indx,indy,group0)
    
    
end

do_pca = 0;
if do_pca
   [coeff,score,latent] = princomp([x(:),y(:)]);
    ax1 = coeff(:,1);% *sqrt(latent(1));
    ax2 = coeff(:,2);% *sqrt(latent(2));
    hold on; plot([0 coeff(1,1)],[0 coeff(2,1)],'r','LineWidth',2);
    hold on; plot([0 coeff(1,2)],[0 coeff(2,2)],'r','LineWidth',2); 
end


% %% Slope test
% if plot_Sch_sens
%     pls_stats2 = pls_stats;
%     pls_stats2(:,group(3).ctgs) = pls_stats(:,group(3).ctgs) ./ pls_stats(:,group(1).ctgs);
%     figure; % subplot(221); 
%     plott_fit(pls_stats2(~bad_any,group(1).ctgs),pls_stats2(~bad_any,group(3).ctgs),'k.');
%     hold on; plot_scattergroups(pls_stats2(~bad_any,group(1).ctgs),pls_stats2(~bad_any,group(3).ctgs),[group(1),group(3)],~bad_any);
%     xlabel('FFC pls(:,1)'); ylabel('FFC pls(:,2)');
% end



end