


function [h1,hstats,pstats,out] = plot_phistats_custstruct(group,opts)
    %% [h1,hstats,pstats] = plot_stats_custstruct(group)
    
    % Set up options struct
    if nargin < 2
        opts = [];
    end
    
    if ~isfield(opts,'paperfig_mode')       % paperfig_mode: 1=clean plots for paper; 0=plots with additional information for viewing
        opts.paperfig_mode = 0;
    end
    
    if ~isfield(opts,'hmask')       % Mask to black out comparisons in statistical testing
        opts.hmask = [];            % A NxN logical matrix, where N=length(group)
    end                             % Set entries to 0 that you don't want to compare
    
    paperfig_mode = opts.paperfig_mode;
    hmask = opts.hmask;

    % Nans are due to their being only 1 data point. Remove these by
    % setting to zero (optional)
    %fign;
    clf;
    
    remove_nans_in_stats = 1;
    use_nonparametric_stats = 1;
        force_unpaired_statistics = 0;
    do_boxplot = 1;
        delete_outliers=0;
        use_numbers_label=1;    % Use numbers to label x-axis of boxplot instead of the legendarr
    do_bonferroni = 1;
    
    
     % Extract data and parametric measures from group
    Ngroups = length(group);
    Xmu = zeros(1,Ngroups);
    Xste = zeros(1,Ngroups);
    h = false(Ngroups,Ngroups);
    p = zeros(Ngroups,Ngroups);
    for i = 1:Ngroups
       X{i} = group(i).datastats(:);
       if ~isreal(X{i}) X{i} = angle(X{i}); end
       Xgroups{i} = i * ones(length(X{i}),1);
       Xmu(i) = circ_mean(X{i});
       Xste(i) = circ_std(X{i});
    end
    Xall = vertcat(X{:});
    Xall_groups = vertcat(Xgroups{:});
    
    
    num_warnings = 0;
    
    statfun = @circ_ktest;
    
    % Calculate p values
    for i = 1:Ngroups
        for j = 1:Ngroups
            [p(i,j)] = statfun(X{i},X{j});
        end
    end
    h = p < 0.05;
    
    
    if ~isempty(hmask)
        if any(size(hmask) ~= size(p)); error('size(hmask) must be NxN, where N is length(group)');end
        p(hmask==0) = 1;    % Failure
        h(hmask==0) = 0;    % False
    end
    
    
    % Do bonferonni if necessary
    hthresh = 0.05;
    Ncomparisons = 1;
    if do_bonferroni;
        Ncomparisons = (Ngroups.^2-Ngroups)/2;  % Number of comparisons is total # elements, minus diagonal, divided by 2
        if ~isempty(hmask)
            htriu = triu(hmask,1);  % Upper triangle of hmask
            Ncomparisons = sum(htriu(:));   % This is the actual # of comparisons
        end
        hthresh = hthresh / Ncomparisons;           % Scale hthresh by appropriate value
        h = p < hthresh;                        % Refresh h;
    end
    
    % Clean up nans
    if remove_nans_in_stats
        nans = isnan(p(:));
        if sum(nans) > 1; warning('Nans present in stats')
        keyboard; end
    
        p(nans) = 1;
        h(nans) = 0;
        
        nans = isnan(Xmu);
        Xmu(nans) = 0;
        Xste(nans) = 0;
    end
    
    panova = anova1(Xall,Xall_groups,'off');
    legend_arr = get_legendarr(group,3,0);
    
    % Plot the bar graph
    if ~do_boxplot
        %[h1, h2] = bar_errsig(h(1,:),Xste,Xmu,'myfontsize',get_fontlevel(1)*get_scaling_factor);    % Bar graph with sig stars (old)
        [h1] = barwitherr(Xste,Xmu,'k');                                                            % Normal bar graph
        %[h1] = bar(Xmu,'k');
    else
        for i = 1:length(X)
            if use_numbers_label; G{i} = i*ones(size(X{i}));
            else G{i} = repmat(legend_arr{i},size(X{i})); end
        end
        Xg = cellfun(@(x) x(:), X,'UniformOutput',false);
        if use_numbers_label; Gg = cellfun(@(x) x(:), G,'UniformOutput',false);
        else Gg = G; end
        [h1] = boxplot( vertcat(Xg{:}), vertcat(Gg{:}),'boxstyle','outline','plotstyle','traditional','extrememode','clip','colors','k');
        if delete_outliers; delete(findobj(gca,'tag','Outliers')); end
        
    end
    
    
    % Add xtick labels (only for bar graph; boxplot turns xtick labels off
    % and uses its own labels. Can't get horizontal ones working for boxplot)
    if  ~do_boxplot
        if any(cellfun(@length,legend_arr) > 3); xticklabel_rotate(1:length(legend_arr),90,legend_arr, 'FontSize',get_fontlevel(1)*get_scaling_factor*0.5)  % If long legends, rotate
        else set(gca,'XTickLabel',legend_arr);                                                      % Else plot normally
        end
    end
    
    % Add some labels
    if ~paperfig_mode
        title(['Anova p=' num2str(panova)]);
        ylabel([group(1).data_name ' phi']);
    end
    fbr=group(1).freqband_stats;
    title(['@ ' num2str(mean(fbr)) ' Hz']);
    set(gca,'Box','off')
    
%     % Enforce test plot
%     clf;
%     bar([1,2,3,4])
%     p=ones(4)*0.04
%     p(1,2) = 0.01; p(1,4) = 0.01;
%     h = p < 0.05;
     
    % Add significance stars
    htop = triu(h,1);                                                 % This chunk of code converts h matrix into 2-element vectors that ID significantly different groups
    [i,j] = ind2sub([Ngroups,Ngroups],find(htop(:)));                 % Get indices of significantly different groups
    matij = [i,j];
    groupij = mat2cell(matij,[ones(1,size(matij,1))],size(matij,2)); % Convert rows of matij into a cell array, one cell for each row
    pij = p(htop(:));
    pij = pij * Ncomparisons;                                   % Hack for bonferroni correction. sigstar forces hthresh to be 0.05; increasing pij is equivalent to decreasing hthresh.
    myhandle = sigstar(groupij,pij);
    
    
    % Rescale axes above stars, as needed
    if ~isempty(groupij)
        if length(groupij) > 1
            yline_max = cellfun(@max,get(myhandle(:,1),'YData'));         % Y location of maximum "line"
            ysigstar = cellfun(@(x) x(2), get(myhandle(:,2),'Position')); % Pull out y position of each "star"
        else
            yline_max = max(get(myhandle(:,1),'YData'));         % Y location of maximum "line"
            ysigstar = get(myhandle(:,2),'Position');            % Pull out y position of each "star"
            ysigstar = ysigstar(2);
        end
        ymax = max([yline_max(:); ysigstar(:)]);	% Get overall max of significance bars
        ymax_curr = ylim;                  % Get current ymax
        if ymax_curr(2) < ymax                % If it's less than sig star level, increase
            ylim([ymax_curr(1) ymax*1.05])
        end
    end
    
    if ~isempty(group(1).ylims_desired); ylim(group(1).ylims_desired); end
    
    hstats = h;
    pstats = p;

    out = [];
end




%     % Test around some syntheitc data
%   
%     mu=1; kk=2.5;
%     y=circ_vmrnd(mu,kk,200);
%     est_var=circ_var(y)         % Empirical variance
%     anal_var=1-besseli(1,kk)/besseli(0,kk) % theoretical variance
%     
%     mu=0; kk=2.5;
%     z=circ_vmrnd(mu,kk,200);
%     circ_kuipertest(y,z,100,1)
%     
%     
%     
%     mu=1; kk=3.5;
%     y=circ_vmrnd(mu,kk,200);
%     mu=0; kk=2.5;
%     z=circ_vmrnd(mu,kk,200);
%     [p] = circ_ktest(y,z)
% %     
    
