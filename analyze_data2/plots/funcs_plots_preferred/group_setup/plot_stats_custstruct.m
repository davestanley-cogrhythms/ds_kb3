

function [handle,hstats,pstats,out,starsum] = plot_stats_custstruct(group,opts)
    %% [handle,hstats,pstats] = plot_stats_custstruct(group)
    
    hstats = [];
    pstats = [];
    
    % Set up options struct
    if nargin < 2
        opts = [];
    end
    
    if ~isfield(opts,'paperfig_mode')       % paperfig_mode: 1=clean plots for paper; 0=plots with additional information for viewing
        opts.paperfig_mode = 0;
    end
    
    if ~isfield(opts,'show_legend')       % paperfig_mode: 1=clean plots for paper; 0=plots with additional information for viewing
        opts.show_legend = false;
    end
    
    if ~isfield(opts,'hmask')       % Mask to black out comparisons in statistical testing
        opts.hmask = [];            % A NxN logical matrix, where N=length(group)
    end                             % Set entries to 0 that you don't want to compare
    
    if ~isfield(opts,'remove_dependent')       % Remove dependent electrode pairs for statistical testing.
        opts.remove_dependent = 0;
    end
    
    if ~isfield(opts,'do_stats')      
        opts.do_stats = 1;
    end
    
    if ~isfield(opts,'do_diagtest')      
        opts.do_diagtest = 0;
    end
    
    if ~isfield(opts,'multicomparisons_mode')      
        opts.multicomparisons_mode = 2;  % 0-None; 1-Bonferroni; 2-Holm-Bonferroni
    end
    
    if ~isfield(opts,'opts_diagtest')      
        opts.opts_diagtest.testmode = 1;        % 1-Binomial test against expected
                                                % 2-Signrank test against threshold
        opts.opts_diagtest.threshold = 1.0;       
    end
    
    % Import data from struct
    paperfig_mode = opts.paperfig_mode;
    hmask = opts.hmask;
    hmask = logical(hmask);
    remove_dependent = opts.remove_dependent;
    show_legend = opts.show_legend;
    do_stats = opts.do_stats;
    do_diagtest = opts.do_diagtest;
        opts_diagtest.testmode = opts.opts_diagtest.testmode;
        opts_diagtest.threshold = opts.opts_diagtest.threshold;
    multicomparisons_mode = opts.multicomparisons_mode;

    % Nans are due to their being only 1 data point. Remove these by
    % setting to zero (optional)
    %fign;
%     clf;
    
    % Define some parameters
    remove_nans_in_stats = 1;
    use_nonparametric_stats = 1;
        force_unpaired_statistics = 0;
    do_boxplot = 0;
        delete_outliers=0;
        use_numbers_label=1;    % Use numbers to label x-axis of boxplot instead of the legendarr
    if do_diagtest
        bino_alpha = 0.05;
    end

    % Set up variables
    Ngroups = length(group);
    
    % Set up hmask
    if isempty(hmask)
        % % % % For neighboring % % % % 
        %hmask = true(size(p)); % For all to all
        % % % % For none to none % % % % 
        %hmask = false(size(p)); % For all to all
        
        hmask = false(Ngroups,Ngroups);
        
        % % % % For make all neighbouring true for normal stats % % % % 
        if do_stats
            myblkdiag = repmat({true(2)},1,floor(Ngroups/2));
            hmask_blk = blkdiag(myblkdiag{:});
            sz = size(hmask_blk);
            hmask(1:sz(1),1:sz(2)) = logical(hmask_blk);
            hmask(logical(eye(size(hmask)))) = false;
        end
        
        % Make diagonal all true for binotest
        if do_diagtest
            hmask(logical(eye(size(hmask)))) = true;
        end
    end
    
    
    if ~isempty(hmask) && (do_stats || do_diagtest)
        if any(size(hmask) ~= [Ngroups,Ngroups]); error('size(hmask) must be NxN, where N is length(group)');end
    end
    
    
    % Extract data and parametric measures from group
    Xmu = zeros(1,Ngroups);
    Xste = zeros(1,Ngroups);
    p = ones(Ngroups,Ngroups);
    for i = 1:Ngroups
       X{i} = group(i).datastats(:);
       Xgroups{i} = i * ones(length(X{i}),1);
       Xmu(i) = mean(X{i});
       if remove_dependent
           Xste(i) = calc_ste_remove_dependent(X{i},group(i).metadata.mypairs);
       else
           Xste(i) = std(X{i}) / sqrt(length(X{i}));
       end
    end
    Xall = vertcat(X{:});
    Xall_groups = vertcat(Xgroups{:});
    
    if do_stats
        [p] = calculate_p_allbars(X,group,remove_dependent,use_nonparametric_stats,force_unpaired_statistics);
        p
    end
    
    if do_diagtest
        switch opts_diagtest.testmode
            case 1
                mytest = @(x) myBinomTest2(x,bino_alpha);
            case 2
                mytest = @(x) signrank(x-opts_diagtest.threshold);
        end
        
        [p_bino] = calculate_p_diagtest(X,group,mytest,remove_dependent);
        p_bino
        ind = logical(eye(size(p)));
        p(ind) = p_bino;
    end
    
    % Apply hmask
    p(hmask==0) = 1;     % Failure
    
    alpha = 0.05;
    [p_bonferroni_sigstars,h1,h2,h3] = bonferroni_for_sigstar(p,alpha,hmask,multicomparisons_mode);
        
    
    
    % Clean up nans
    if remove_nans_in_stats
        nans = isnan(p(:));
        if sum(nans) > 1; warning('Nans present in stats')
        keyboard; end
    
        p(nans) = 1;
        p_bonferroni_sigstars(nans) = 1;
        
        nans = isnan(Xmu);
        Xmu(nans) = 0;
        Xste(nans) = 0;
    end
    
    panova = anova1(Xall,Xall_groups,'off');
    legend_arr = get_legendarr(group,3,0);
    
    % Plot the bar graph
    if ~do_boxplot
        %[h1, h2] = bar_errsig(h(1,:),Xste,Xmu,'myfontsize',get_fontlevel(1)*get_scaling_factor);    % Bar graph with sig stars (old)
        [handle] = barwitherr(Xste,Xmu,'k');                                                            % Normal bar graph

        
    else
        for i = 1:length(X)
            if use_numbers_label; G{i} = i*ones(size(X{i}));
            else G{i} = repmat(legend_arr{i},size(X{i})); end
        end
        Xg = cellfun(@(x) x(:), X,'UniformOutput',false);
        if use_numbers_label; Gg = cellfun(@(x) x(:), G,'UniformOutput',false);
        else Gg = G; end
        [handle] = boxplot( vertcat(Xg{:}), vertcat(Gg{:}),'boxstyle','outline','plotstyle','traditional','extrememode','clip','colors','k');
        if delete_outliers; delete(findobj(gca,'tag','Outliers')); end
        
    end
    
    
    % Add xtick labels (only for bar graph; boxplot turns xtick labels off
    % and uses its own labels. Can't get horizontal ones working for boxplot)
    if  ~do_boxplot
        if any(cellfun(@length,legend_arr) > 3);

            if is_newgraphics_matlab
                set(gca,'XTickLabel',legend_arr,'XTickLabelRotation',90,'TickLabelInterpreter',group(1).interpreter);
            else
                xticklabel_rotate(1:length(legend_arr),90,legend_arr, 'FontSize',get_fontlevel(1)*get_scaling_factor*0.5,'interpreter',group(1).interpreter)  % If long legends, rotate - without interpreter
            end
            
        else 
            temp=legend_arr;
            temp=strrep(temp,'_','\_');
            set(gca,'XTickLabel',temp);                                                     % Else plot normally
        end
    end
    
    % Add some labels
    if ~paperfig_mode
        title(['Anova p=' num2str(panova)]);
        ylabel(group(1).data_name);
    end
    
    
    fbr=group(1).freqband_stats;
    if ~isempty(fbr) && ~strcmp(group(1).data_type,'FR');
        title(['@ ' num2str(mean(fbr)) ' Hz']);
    end
    
    % Adjust position of axis to make sure everything fits
    
%     if ~is_newgraphics_matlab
%         temp=get(gca,'Position');
%         temp(4)=temp(4)-0.05;
%         set(gca,'Position',temp);
%     end
    set(gca,'Box','off')
    
    
%     % Enforce test plot
%     clf;
%     bar([1,2,3,4])
%     p=ones(4)*0.04
%     p(1,2) = 0.01; p(1,4) = 0.01;
%     h = p < 0.05;
     
    % % % % Add significance stars % % % % 
        % This chunk of code converts h matrix into 2-element vectors that
        % ID significantly different groups (groupij) and their
        % corresponding p values.
    htop = triu(p_bonferroni_sigstars < alpha,0);                     % P values with bonferroni correction.                
    [i,j] = ind2sub([Ngroups,Ngroups],find(htop(:)));                 % Get indices of significantly different groups
    matij = [i,j];
    groupij = mat2cell(matij,[ones(1,size(matij,1))],size(matij,2)); % Convert rows of matij into a cell array, one cell for each row
    pij = p_bonferroni_sigstars(htop(:));
        % Now do sig stars
    myhandle = sigstar(groupij,pij);
    hstats = p_bonferroni_sigstars < alpha;
    pstats = p;
    % % % % End significance stars % % % % 
    
    
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
        ylims_curr = ylim;                  % Get current ymax
        if ylims_curr(2) < ymax                % If it's less than sig star level, increase
            ylim([ylims_curr(1) ymax*1.05])
        end
    end
    
    if ~isempty(group(1).zlims_desired); ylim(group(1).zlims_desired); end
    if show_legend;
        try
            [numcells] = arrayfun(@(s) sum(s.numcells),group);
        catch
            numcells = group.myfunc_out(@(s) sum(s.numcells));
        end
        if length(unique(numcells)) == 1
            legend(['N=' num2str(numcells(1))]);
        end
    end
    
    
    % Add horizontal line for number expected by chance
    if do_diagtest
        numcells_all = [group.numcells];
        
 
        
        switch opts.opts_diagtest.testmode
            case 1; yval = [bino_alpha*numcells_all(1)];
            case 2; yval = opts.opts_diagtest.threshold;
        end
        
        if (opts.opts_diagtest.testmode == 1 && all(numcells_all == numcells_all(1))) ...   % This line only makes sense if # expected is same for all bars
                || (opts.opts_diagtest.testmode == 2 && yval > 0)             
            xlims_curr = xlim;
            hold on; plot([xlims_curr(1), xlims_curr(2)],[yval, yval],'--','LineWidth',3,'Color',[0.5,0.5,0.5]);
        end
    end
    

%     index = f >= fb(1) & f < fb(2);
%     X1=(group1);
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


function Xout = calc_ste_remove_dependent(Xin,mypairsx)
    
    bads = ~remove_dependent_FFC_pairs(mypairsx);
    Xout = std(Xin(~bads)) / sqrt(length(Xin(~bads)));
end

function [p,h] = stats_remove_dependent_paired(X,Y,mypairsx,mypairsy,mystatfunc)
    bads = ~remove_dependent_FFC_pairs(mypairsx);
    
    [p,h] = mystatfunc(X(~bads),Y(~bads));
    
    
end

function [p,h] = stats_remove_dependent_unpaired(X,Y,mypairsx,mypairsy,mystatfunc)
    bads = ~remove_dependent_FFC_pairs([mypairsx; mypairsy]);
    badsx = bads(1:length(X));
    badsy = bads(length(X)+1:(length(X)+length(Y)));
    
    X = X(~badsx);
    Y = Y(~badsy);
    
    [p,h] = mystatfunc(X,Y);
    
end


function ispaired = cells_are_paired(cells1,cells2)
    ispaired = 0;
    if size(cells1) == size(cells2)
        ispaired = all(cells1(:) == cells2(:));
    end
end


function [p,h1] = calculate_p_allbars(X,group,remove_dependent,use_nonparametric_stats,force_unpaired_statistics)

    % Extract data and parametric measures from group
    Ngroups = length(group);
    h1 = false(Ngroups,Ngroups);
    p = ones(Ngroups,Ngroups);
    

    num_warnings = 0;
    
    % Calculate p values
    for i = 1:Ngroups
        for j = 1:Ngroups
            if ~use_nonparametric_stats
                if cells_are_paired(group(i).cells(:),group(j).cells(:)) && ~force_unpaired_statistics   % DETECTED PAIRED
                        [h1(i,j), p(i,j)] = ttest(X{i},X{j}); % Assume data paired
                        if remove_dependent; [h1(i,j), p(i,j)] = stats_remove_dependent_paired(X{i},X{j},group(i).metadata.mypairs,group(j).metadata.mypairs,@ttest); end
                else
                    [h1(i,j), p(i,j)] = ttest2(X{i},X{j});
                    if remove_dependent; [h1(i,j), p(i,j)] = stats_remove_dependent_unpaired(X{i},X{j},group(i).metadata.mypairs,group(j).metadata.mypairs,@ttest2); end
                    if force_unpaired_statistics && num_warnings == 0; warning('Using unpaired when really should use paired'); num_warnings=num_warnings+1; end
                end
            else
                if cells_are_paired(group(i).cells(:),group(j).cells(:)) && ~force_unpaired_statistics   % DETECTED PAIRED
                        [p(i,j), h1(i,j)] = signrank(X{i},X{j}); % Assume data paired
                        if remove_dependent; [p(i,j), h1(i,j)] = stats_remove_dependent_paired(X{i},X{j},group(i).metadata.mypairs,group(j).metadata.mypairs,@signrank); end
                else
                    [p(i,j), h1(i,j)] = ranksum(X{i},X{j});
                    if remove_dependent; [p(i,j), h1(i,j)] = stats_remove_dependent_unpaired(X{i},X{j},group(i).metadata.mypairs,group(j).metadata.mypairs,@ranksum); end
                    if force_unpaired_statistics && num_warnings == 0; warning('Using unpaired when really should use paired'); num_warnings=num_warnings+1; end
                end
            end
        end
    end
    
    
end

function [p_bino] = calculate_p_diagtest(X,group,mytest,remove_dependent)
    % bino_alpha = 0.05;
    
    Ngroups = length(group);
    
    p_bino = zeros(1,Ngroups);
    for i = 1:Ngroups
        
        Xcurr = X{i};
        if remove_dependent
            bads = ~remove_dependent_FFC_pairs(group(i).metadata.mypairs);
            Xcurr = Xcurr(~bads);
        end
        
        p_bino(i) = mytest(Xcurr);
        
        %[s/n p_bino(i)]

    end
    
end

function p_bino = myBinomTest2(Xcurr,bino_alpha)

        s = sum(Xcurr > 0.5);
        n = length(Xcurr);
        p = bino_alpha;
        sided = 'one';
        [p_bino] = myBinomTest(s,n,p,sided);

end

function [pout,h1,h2,h3] = bonferroni_for_sigstar(p,alpha,hmask,multicomparisons_mode)

    % Do bonferonni if necessary
    Ngroups = length(p);
    
    %Ncomparisons = (Ngroups.^2-Ngroups)/2;  % Number of comparisons is total # elements, minus diagonal, divided by 2
    Ncomparisons = (Ngroups.^2-Ngroups)/2 + Ngroups;  % We're now allowing comparisons along the diagnoal. 
    
    % Set up pout
    pout = p;
    
    if ~isempty(hmask)    
        htriu = triu(hmask,0);  % Upper triangle of hmask
        Ncomparisons = sum(htriu(:));   % This is the actual # of comparisons
        if Ncomparisons < 1; Ncomparisons = 1; end
    end
    if multicomparisons_mode == 0
        pout = p;
        h1 = p < alpha;
        h2 = p < alpha / 5;
        h3 = p < alpha / 50;
    elseif multicomparisons_mode == 1
        h1 = p < alpha/Ncomparisons;                       % Refresh h;
        h2 = p < alpha/5/Ncomparisons;                     % Refresh h;
        h3 = p < alpha/50/Ncomparisons;                    % Refresh h;
        if sum(h1) < sum(h2); warning('h1 missing rejections'); end
        if sum(h2) < sum(h3); warning('h2 missing rejections'); end
        pout = ones(size(p));
        pout(h1) = 0.04;   % Force values of p for plotting sigstar!
        pout(h2) = 0.009;
        pout(h3) = 0.0009;
    elseif multicomparisons_mode == 2
        [~, p1] = holm_bonferroni(p(htriu),alpha);  % Get critical p value for chosen data.
        [~, p2] = holm_bonferroni(p(htriu),alpha/5);
        [~, p3] = holm_bonferroni(p(htriu),alpha/50);
        if p2 > p1; warning('p2 > p1'); end
        if p3 > p2; warning('p3 > p2'); end
        h1 = p <= p1;       % Need to do this to restore shape
        h2 = p <= p2;
        h3 = p <= p3;
        if sum(h1) < sum(h2); warning('h1 missing rejections'); end
        if sum(h2) < sum(h3); warning('h2 missing rejections'); end
        pout = ones(size(p));
        pout(h1) = 0.04;   % Force values of p for plotting sigstar!
        pout(h2) = 0.009;
        pout(h3) = 0.0009;
    end
    

end


