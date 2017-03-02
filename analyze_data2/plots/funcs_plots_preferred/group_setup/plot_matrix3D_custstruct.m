


function [h1, out] = plot_matrix3D_custstruct(abscissa,group,opts_PM3D,opts,my_clist,my_linestyle,my_linewidth)

% Set defaults
    if nargin < 7
        my_linewidth = [];
    end
    if nargin < 6
        my_linestyle = [];
    end
    
    if nargin < 5;
        my_clist = [];
    end
    
    if nargin < 4
        opts = [];
    end
    
    if nargin < 3
        opts_PM3D = [];
    end
    
    % Fill in empty fields
    if isempty(opts_PM3D); opts_PM3D = {'do_mean',1,'do_zscore',0,'showErrorbars',1}; end
    
    if ~isfield(opts,'paperfig_mode'); opts.paperfig_mode = 0; end % paperfig_mode: 1=clean plots for paper; 0=plots with additional information for viewing
    if ~isfield(opts,'do_sgolay'); opts.do_sgolay = 0; end
    if ~isfield(opts,'stats_mode'); opts.stats_mode = 0; end        % stats_mode: 2=different from 1; 3=different from 0
    if ~isfield(opts,'remove_dependent'); opts.remove_dependent = 0; end
    if ~isfield(opts,'do_subplots'); opts.do_subplots = 0; end
    if ~isfield(opts,'groups_per_subplot'); opts.groups_per_subplot = 4; end
    if ~isfield(opts,'max_subplots_per_fig'); opts.max_subplots_per_fig = 16; end
    
    if isempty(my_clist); my_clist = get_clist; end
    if isempty(my_linestyle); my_linestyle = repmat({'-'},1,length(group)); end
    if isempty(my_linewidth); my_linewidth = repmat(2,1,length(group)); end
    
    do_subplots = opts.do_subplots;
        max_subplots_per_fig = opts.max_subplots_per_fig;
        use_subplot_grid = 1;
    groups_per_subplot = opts.groups_per_subplot;
        
        
    if do_subplots
        curr_subplots = 1;
        hsp = [];
        N_subplots = ceil(length(group)/groups_per_subplot);
        max_subplots_per_fig = min(max_subplots_per_fig,N_subplots);         % Reduce max subplots if length of group is less
        for i = 1:N_subplots
            % Manage figure creation
            [hsp, curr_subplots,returns] = new_subplot(max_subplots_per_fig,curr_subplots,hsp,use_subplot_grid);    % Create a new subplot entry or start a new figure if current fig is full.
            if strcmp(returns,'q') || strcmp(returns,'Q'); break; end
            
            % Import some metadata info from group(1)
            group_curr = group((i-1)*groups_per_subplot+1:i*groups_per_subplot);
            group_curr(1).xlims_desired = group(1).xlims_desired;
            group_curr(1).zlims_desired = group(1).zlims_desired;
            group_curr(1).data_name = group(1).data_name;
            [h1, out] = plot_matrix3D_custstruct_subfunc(abscissa,group_curr,opts_PM3D,opts,my_clist,my_linestyle,my_linewidth);
            set(gca,'FontSize',10);
        end
        
    else
        [h1, out] = plot_matrix3D_custstruct_subfunc(abscissa,group,opts_PM3D,opts,my_clist,my_linestyle,my_linewidth);
    end


end



function [h1, out] = plot_matrix3D_custstruct_subfunc(abscissa,group,opts_PM3D,opts,my_clist,my_linestyle,my_linewidth)


    % Import variables from structures
    paperfig_mode = opts.paperfig_mode;
    do_sgolay = opts.do_sgolay;
    stats_mode = opts.stats_mode;
    remove_dependent = opts.remove_dependent;
        mean_frequencies = 0;       % 0 - do stats on a frequency range
                                    % 1 - do stats on a single point in frequency space
    
    
                                    
    % First loop to remove nans and infs from data
    for i = 1:length(group);
        nanind = isnan(group(i).data);
        found_bads = 0;
        if any(nanind(:))
            found_bads = 1;
            group(i).data(nanind) = 1;
            warning('NaNs detected in group. Converting to ones.');
        end
        infind = isinf(group(i).data);
        if any(infind(:))
            found_bads = 1;
            group(i).data(infind) = 1;
            warning('Infs detected in group. Converting to ones.');
        end
        
        if found_bads
            group(i).data = group(i).data(2:end-1,:);
            group(i).xdata = group(i).xdata(2:end-1);
            abscissa = abscissa(2:end-1);
        end
        
    end
    
    
    % Calculate means across cells if not already populated
    for i = 1:length(group)
        
        if remove_dependent
            bad_dependent = ~remove_dependent_FFC_pairs(group(i).metadata.mypairs);
        else
            bad_dependent = false(1,size(group(i).data,2));
        end
        
        % Calculate mean of data
        if isempty(group(i).data_mu)
            group(i).data_mu = mean(group(i).data,2);
        end
        
        % Calculate STE of data
        if isempty(group(i).data_STE)
            group(i).data_STE=std(group(i).data(:,~bad_dependent),[],2) / sqrt(size(group(i).data(:,~bad_dependent),2)); % Standard error  - only across non-dependent
        end
        
    end

    

    
    % Plot the data
    h1=[];
    for i = 1:length(group)
        [data_mu, x] = extract_data(group(i),abscissa,do_sgolay,[],[],'data_mu');
        
        if ~isempty(group(i).colour_desired)
            colorspec = {'Color',group(i).colour_desired};
        else
            if ischar(my_clist); colorspec = {'Color',my_clist(i)};    % If color list is in char
            else colorspec = {'Color',my_clist(i,:)};                % If color list is a vector
            end
        end
        
        if ~isempty(data_mu)            
            %hold on; [h1(i)] = plott_matrix3D(x,data,opts_PM3D{:},'LineSpec',{colorspec{:},'LineWidth',2});   % Old code, using plott_matrix3D. This code is okay in most cases,
                                                                                                                % but doesn't work with remove_dependent turned on.
                                                                                                                
                                                                                                                
                                                                                                                
            data_STE = group(i).data_STE;
            LineSpec1 = {colorspec{:},'LineWidth',my_linewidth(i),'LineStyle',my_linestyle{i}};
            if ischar(colorspec{2});
                LineSpec2 = [colorspec{2}, my_linestyle{i}];
                cmap = [];
            else
                LineSpec2 = [my_linestyle{i}];
                cmap = colorspec{2};
            end
            hold on; [h1(i)] = plott_matrix3D_simplified(x,data_mu,data_STE,LineSpec1,LineSpec2,cmap,opts_PM3D{find(strcmp(opts_PM3D,'showErrorbars'))+1});
            
                                                                                                                % ^^ This is a really cheap temporary hack!
            
        else
            fprintf('Empty group \n');
            hold on; [h1(i)] = plot(any(nanind,2),colorspec{:});       % Plot a dummy line to allow for legend entry
        end
        
    end
    
    legend_arr = get_legendarr(group);
    
   
    
    % % % % Do statistical testing
    if ~paperfig_mode && size(data,2) > 1
        if length(group) > 2            % If > 2 groups, do anova
            Ngroups = length(group);
            Xmu = zeros(1,Ngroups);
            Xste = zeros(1,Ngroups);
            h = zeros(Ngroups,Ngroups);
            p = h;
            for i = 1:Ngroups
               X{i} = group(i).datastats(:);
               Xgroups{i} = i * ones(length(X{i}),1);
               Xmu(i) = mean(X{i});
               Xste(i) = std(X{i}) / sqrt(length(X{i}));
            end
            Xall = vertcat(X{:});
            Xall_groups = vertcat(Xgroups{:});


            if ~isreal(Xall)
                Xall = pls_complex2angle_centpi(Xall);
                warning('Complex data received! We should be bootstrapping now, meaning that input data should always be real.');
            end

            panova = anova1(Xall,Xall_groups,'off');
            title(['Anova p=' num2str(panova) ' for freqband ' num2str(group(1).freqband_stats(1)) '-' num2str(group(1).freqband_stats(2)) 'Hz']); 
        elseif length(group) > 1   % If only 2 groups, do ttest / ranksum
    %         [h,p] = ttest(group(1).datastats(:),group(2).datastats(:));title(['t test p=' num2str(p) ' for freqband ' num2str(group(1).freqband_stats(1)) '-' num2str(group(1).freqband_stats(2)) 'Hz']); 
    %         [h,p] = ttest2(group(1).datastats(:),group(2).datastats(:));title(['t test p=' num2str(p) ' for freqband ' num2str(group(1).freqband_stats(1)) '-' num2str(group(1).freqband_stats(2)) 'Hz']); 
            if all(size(group(1).datastats(:)) == size(group(2).datastats(:)))
                [p,h] = signrank(group(1).datastats(:),group(2).datastats(:)); title(['Signrank p=' num2str(p) ' for freqband ' num2str(group(1).freqband_stats(1)) '-' num2str(group(1).freqband_stats(2)) 'Hz']); % PAIRED
            else
                [p,h] = ranksum(group(1).datastats(:),group(2).datastats(:)); title(['Ranksum p=' num2str(p) ' for freqband ' num2str(group(1).freqband_stats(1)) '-' num2str(group(1).freqband_stats(2)) 'Hz']); % UNPAIRED       
            end
        end
    end
    
    
    if isempty(strfind(group(1).xdata_name,'Time'))
        if isempty(group(1).zlims_desired)
            [myylims] = ylim;
        else
            myylims = group(1).zlims_desired;
        end
        hold on; plot([11,11],[myylims],'k--','LineWidth',2)
        hold on; plot([18,18],[myylims],'k-','LineWidth',2)
    else
        add_stages_separators([-1 0 2 3 6 ],'k'); add_stages_names([-1 0 2 3 6 ],'k'); 
    end
%     hold on; plot([26 26],[myylims],'k-','LineWidth',2)
    
    % Arrange plot
    if ~isempty(h1)
        if is_newgraphics_matlab
            legend([h1],legend_arr,'Interpreter',group(1).interpreter,'Location','NorthEast');
        else
            legend([h1],legend_arr,'Interpreter',group(1).interpreter,'Location','NorthEast')
        end
    end
    
    if ~isempty(group(1).data_name)
        if ~paperfig_mode
            ylabel(group(1).data_name(1,:));
        end
    end

    % Adjust position of axis to make sure everything fits
        temp=get(gca,'Position');
        temp(2)=temp(2)+0.05; temp(4)=temp(4)-0.05;
        set(gca,'Position',temp);

    xlabel(group(1).xdata_name(1,:));
    
    if ~isempty(group(1).zlims_desired); ylim(group(1).zlims_desired); end
    if ~isempty(group(1).xlims_desired); xlim(group(1).xlims_desired); end
    
%     xlim_temp = xlim;
%     if xlim_temp <= 100
%         set(gca,'XTick',0:10:xlim_temp(2));     % Optional - enforce tic spacing by 10
%     elseif xlim_temp <= 150
%         set(gca,'XTick',0:10:xlim_temp(2));     % Optional - enforce tic spacing by 10
%     elseif xlim_temp <= 200
%         set(gca,'XTick',0:20:xlim_temp(2));     % Optional - enforce tic spacing by 10
%     elseif xlim_temp <= 20
%         set(gca,'XTick',0:25:xlim_temp(2));     % Optional - enforce tic spacing by 10
%     end

    
    
    % Plotting of significance dots!
    switch stats_mode
        case 0   % Do nothing
        case 1 
            warning('Not implemented');
        case {2,3,4}
            switch stats_mode
                case 2
                    different_from = 1;
                    star_spacing = 0.02;
                    group_spacing =1;
                case 3
                    different_from = 0;
                    star_spacing = 0.002;
                    group_spacing = 1;
                case 4
                    group_spacing = 2;
                    star_spacing = 0.02;
            end
            for i = 1:group_spacing:length(group)
                
                [data, x] = extract_data(group(i),abscissa,do_sgolay,[],[],'data');
                if stats_mode == 4
                    [data2, x] = extract_data(group(i+1),abscissa,do_sgolay);
                end
                dx = mode(diff(x));
                if remove_dependent
                    bad_dependent = ~remove_dependent_FFC_pairs(group(i).metadata.mypairs);
                else
                    bad_dependent = false(size(data,2),1);
                end
                
                if ischar(my_clist); colorspec = {'Color',my_clist(i)};    % If color list is in char
                else colorspec = {'Color',my_clist(i,:)};                % If color list is a vector
                end
                
        
                if ~isempty(data)
                    if isempty(strfind(group(1).xdata_name,'Time'))   % FFC or PSD
                        siglimits = [8 100];
                        sample_spacing = round(1.5/dx);           % For tapers mode 4 [TW=3,K=5], bandwidth is 3/0.5=6 or 3/0.8=3.75
                    else                                        % Firing rate
                        siglimits=[min(x) max(x)];
                        sample_spacing = round(0.05/dx);       % For tapers mode=2, should be 0.3. I'm oversampling here for visualization purposes
                    end
                    
                    if ~mean_frequencies
                        inds = find(x > siglimits(1) & x <= siglimits(2));
                        inds = inds(1:sample_spacing:end);
                    else
                        freqband_stats = group(i).freqband_stats;
                        inds = find_closest(x,mean(freqband_stats));
                    end
                    sig = zeros(1,length(inds));
                    for j = 1:length(inds)
                        %sig(j) = signrank(data(inds(j),~bad_dependent)-different_from);
                        if stats_mode ~=4
                            [~,sig(j)] = ttest(data(inds(j),~bad_dependent)-different_from);
                        elseif stats_mode == 4
                            [~,sig(j)] = ttest(data(inds(j),~bad_dependent)-data2(inds(j),~bad_dependent));
                        end
                        
                    end
                    sig2 = sig;
                    sig2 = sig2 * length(inds); % Bonferonni correction
                    thresh = 0.0000001;
                    thresh = 0.05;
                    
                    if ~mean_frequencies
                        myylims = ylim;
                        chosen = inds(find(sig2 < thresh));
                        hold on; plot(x(chosen),(myylims(1))*ones(1,length(chosen))+star_spacing*i,'.','LineWidth',2,'MarkerSize',get_fontlevel(4)*get_scaling_factor,colorspec{:});
                    else
                        if sig2 < thresh; hold on; plot([freqband_stats(1) freqband_stats(2)],1.1*ones(1,2)+0.01*i,'LineWidth',2,colorspec{:}); end
                    end
                end
            end
            
    end
    
%     
%     index = f >= fb(1) & f < fb(2);
%     X1=(group);
%     X2=(group2);
%     if stats_dozscore
%         X1=zscore(X1);
%         X2=zscore(X2);
%     end
%     X1 = mean(X1(index,:,:),1);
%     X1 = squeeze(X1); X1=X1(:);
%     X2 = mean(X2(index,:,:),1);
%     X2 = squeeze(X2); X2=X2(:);
% 
%     
%     Xmu = [mean(X1) mean(X2)];
%     Xste = [ste(X1,[],1) ste(X2,[],1)];
%     [h p] = ttest2(X1,X2)
%         
%     if statsmode == 1
% 
% 
%         [h1 h2 ] = bar_errsig(h,Xste,Xmu);
%         ylabel(['SFC in Range ' num2str(fb(1)) '-' num2str(fb(2)) 'Hz range']);
%     elseif statsmode == 2
%         nbins = 20;
%         subplot(121); hist(X1,nbins);
%         subplot(122); hist(X2,nbins);
%         h1=[];h2=[];
%     end
%     
%     title(['p = ' num2str(p(1))])
%     
    out = [];
end


function [hl1,hl2] = plott_matrix3D_simplified(t,X,Xste,LineSpec,LineSpec2,cmap,show_errorbars)
            t=t(:);
            t = mean(t,2);      % Some extra operations, but it's worth it.
            

            hold on; hl1 = plot(t,squeeze(X),LineSpec{:}); hold on; 
            if show_errorbars
                if isempty(cmap)
                    hold on; hl2 = boundedline(t,squeeze(X),permute(squeeze(Xste),[1 3 2]),LineSpec2,'alpha');
                else
                    hold on; hl2 = boundedline(t,squeeze(X),permute(squeeze(Xste),[1 3 2]),LineSpec2,'alpha',cmap);
                end
            end

end

