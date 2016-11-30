

function [handle,hstats,pstats,out,starsum] = plot_reg_custstruct(rs,opts)
    %% [h1,hstats,pstats] = plot_stats_custstruct(group)
    
    % Set up options struct
    if nargin < 2
        opts = [];
    end
    
    % Plot old code for testing (be sure to clear figure after!)
    clf
    [~,~,~,~,starsum_old] = plot_reg_custstruct_old(rs,opts);    
    clf
    
    if ~isfield(opts,'paperfig_mode')       % paperfig_mode: 1=clean plots for paper; 0=plots with additional information for viewing
        opts.paperfig_mode = 0;
    end
    
    if ~isfield(opts,'hmask')       % Mask to black out comparisons in statistical testing
        opts.hmask = [];            % A NxN logical matrix, where N=length(group)
    end                             % Set entries to 0 that you don't want to compare
    
    paperfig_mode = opts.paperfig_mode;
    hmask = opts.hmask;
    hmask = logical(hmask);

    % Nans are due to their being only 1 data point. Remove these by
    % setting to zero (optional)
    %fign;
    clf;
    
    multicomparisons_mode = 2;  % 0-None; 1-Bonferroni; 2-Holm-Bonferroni
    
    % Extract data and parametric measures from group
    Xmu = [rs.x];
    Xste = [rs.xerr];
    h = [rs.h];
    p = [rs.p];
    legend_arr = {rs.legend_arr};
    
    % Identify things to remove
    ignoreind = strcmp(legend_arr,'(Intercept)') | strcmp(legend_arr,'adj11');
    
    % Remove ignoring entries
    Xmu=Xmu(~ignoreind);
    Xste=Xste(~ignoreind);
    h=h(~ignoreind);
    p=p(~ignoreind);
    legend_arr=legend_arr(~ignoreind);
    
    if isempty(hmask)
        hmask = true(size(p));
    end
    if ~isempty(hmask)
        if any(numel(hmask) ~= numel(p)); error('size(hmask) must be NxN, where N is length(group)');end
        p(hmask==0) = 1;    % Failure
        h(hmask==0) = 0;    % False
    end
    
    % Do bonferonni if necessary
    Ncomparisons = length(p);
    if ~isempty(hmask)
        Ncomparisons = sum(hmask);
    end
    
    alpha = 0.05;
    switch multicomparisons_mode
        case 0
            h2 = zeros(size(p));
            h2(double(p) < alpha) = 1;
            h2(double(p) < alpha/5) = 2;
            h2(double(p) < alpha/50) = 3;
        case 1
            h2 = zeros(size(p));
            h2(double(p) < alpha/Ncomparisons) = 1;
            h2(double(p) < alpha/Ncomparisons/5) = 2;
            h2(double(p) < alpha/Ncomparisons/50) = 3;
        case 2
            [~, p1] = holm_bonferroni(p(hmask),alpha);  % Get critical p value for chosen data.
            [~, p2] = holm_bonferroni(p(hmask),alpha/5);
            [~, p3] = holm_bonferroni(p(hmask),alpha/50);
            if p2 > p1; warning('p2 > p1'); end
            if p3 > p2; warning('p3 > p2'); end
            h2 = zeros(size(p));
            h2(p <= p1)=1;       % Need to do this to restore shape
            h2(p <= p2)=2;
            h2(p <= p3)=3;
    end
    
    % Plot the bar graph
    [handle, herr] = bar_errsig(h2,Xste,Xmu,'k','myfontsize',30,'star_shift_factor',1.1);
    if any(cellfun(@length,legend_arr) > 3); xticklabel_rotate(1:length(legend_arr),90,legend_arr, 'FontSize',get_fontlevel(1)*get_scaling_factor*0.5,'interpreter','none')  % If long legends, rotate
    else
        temp=legend_arr;
        temp=strrep(temp,'_','\_');
        set(gca,'XTickLabel',temp);                                                     % Else plot normally
    end

    
    
    % Add some labels
    if ~paperfig_mode
        ylabel(char(rs(1).Formula));
    end
    ylabel(char(rs(1).Formula));
    fbr=rs(1).freqband_stats;
    
    title(['@ ' num2str(mean(fbr)) ' Hz']);
    set(gca,'Box','off')
    

    
%     % Add significance stars
%     pij = p;
%     pij = pij * Ncomparisons;                                   % Hack for bonferroni correction. sigstar forces hthresh to be 0.05; increasing pij is equivalent to decreasing hthresh.
%     
%     N=length(h);
%     htop = eye(N);                                                 % This chunk of code converts h matrix into 2-element vectors that ID significantly different groups
%     [i,j] = ind2sub([N,N],find(htop(:)));                 % Get indices of significantly different groups
%     matij = [i,j];
%     groupij = mat2cell(matij,[ones(1,size(matij,1))],size(matij,2));
%     
%     warning('Sig star in wrong place');
%     myhandle = sigstar(groupij,pij);
%     
    
%     % Rescale axes above stars, as needed
%     if ~isempty(groupij)
%         if length(groupij) > 1
%             yline_max = cellfun(@max,get(myhandle(:,1),'YData'));         % Y location of maximum "line"
%             ysigstar = cellfun(@(x) x(2), get(myhandle(:,2),'Position')); % Pull out y position of each "star"
%         else
%             yline_max = max(get(myhandle(:,1),'YData'));         % Y location of maximum "line"
%             ysigstar = get(myhandle(:,2),'Position');            % Pull out y position of each "star"
%             ysigstar = ysigstar(2);
%         end
%         ymax = max([yline_max(:); ysigstar(:)]);	% Get overall max of significance bars
%         ymax_curr = ylim;                  % Get current ymax
%         if ymax_curr(2) < ymax                % If it's less than sig star level, increase
%             ylim([0 ymax*1.05])
%         end
%     end
    
    
    hstats = h;
    pstats = p;

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
    
    
    % Check starsum - this slows things down, but it checks to make sure
    % the latest version of code matches up to old version. Once confident
    % it's working, can remove for a speed boost.
    starsum(1) = sum(h2 == 1);
    starsum(2) = sum(h2 == 2);
    starsum(3) = sum(h2 == 3);
    
    if multicomparisons_mode == 1
        if any(starsum_old ~= starsum); warning('Suspected output mismatch between two versions of code'); keyboard; end
    elseif multicomparisons_mode == 2
            % Mode 2 should be at least as good as mode 1. If not,
            % something is wrong
        if any(starsum_old > starsum); warning('Suspected output mismatch between two versions of code'); keyboard; end
    end
    
end