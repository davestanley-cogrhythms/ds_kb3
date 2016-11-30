

function [h1] = plot_histplot(varargin)
    % Megafunction for plotting histograms based on data in group structure
    global bar_orientation


    % Pull out inputs   
    [reg, args]=parseparams(varargin);  % Regs should contain t and X
    
    p = inputParser;
        default_plot_type = 'hist_all';
        valid_plot_typees = {'hist','histall_imagesc','histall_pdf','histall_cdf','hist_paired','hist_autogrouped'};
        check_plot_type = @(x) any(validatestring(x,valid_plot_typees));
    addOptional(p,'plot_type',default_plot_type,check_plot_type);
        default_orient = 'vert';
        valid_orient = {'vert','horiz'};
        check_orient = @(x) any(validatestring(x,valid_orient));
    addOptional(p,'bar_orientation',default_orient,check_orient);
    addOptional(p,'hist_settings',struct,@isstruct);
    
    parse(p,args{:});
    
    % Pull out inputs
    group = reg{1};
    hist_settings = p.Results.hist_settings;
    plot_type = p.Results.plot_type;
    bar_orientation = p.Results.bar_orientation;
    
    % Set up default settings to hist_settings structure
    if ~isfield(hist_settings,'nbins_or_edges'); hist_settings.nbins_or_edges = []; end             % Nbins or edges input to histcounts function
    if ~isfield(hist_settings,'name_value_pairs'); hist_settings.name_value_pairs = {}; end         % Name value pair inputs to histcounts function
    if ~isfield(hist_settings,'show_vertical_bar'); hist_settings.show_vertical_bar = false; end
    if ~isfield(hist_settings,'symmetric_about_x_axis'); hist_settings.symmetric_about_x_axis = false; end
    
    
    % Nans are due to their being only 1 data point. Remove these by
    % setting to zero (optional)
    remove_nans_in_stats = 1;
    
    
    switch plot_type
        case {'hist','histall_imagesc','histall_pdf','histall_cdf'}
            [h1] = hist_all(group,plot_type);
        case 'hist_paired'
            legend_arr = get_legendarr(group);
            N=length(group);
            if mod(N,2) ~= 0
                error('For paired mode, N must be even');
            end

            
            for i = 1:floor(N/2);
                if N>2; figure; end
                [h1] = hist_paired(group(i*2-1).datastats,group(i*2).datastats,hist_settings);
                legend(legend_arr((i*2-1):(i*2)));
            end
            
        case 'hist_autogrouped'
            [h1] = hist_grouped(group,hist_settings);
    end
    
    
    if hist_settings.show_vertical_bar
        yl = ylim;
        hold on; plot([0,0],[yl(1) yl(2)],'k','LineWidth',2);
        
    end
    
    if hist_settings.symmetric_about_x_axis
        xl = xlim;
        xlmax = max(abs(xl));
        xlim([-xlmax xlmax]);
    end
    
    xlabel(group(1).data_name);
    ylabel('N');
end





function [h1] = hist_paired(data1,data2,hist_settings)
    % Plot stacked histogram comparing 2 groups
    
    nbins_or_edges = hist_settings.nbins_or_edges;
    name_value_pairs = hist_settings.name_value_pairs;
    
    i=1; j=2;
    % Estimate best edges for all data
    if ~isempty(nbins_or_edges)
        [~, edges] = histcounts([data1, data2],nbins_or_edges,name_value_pairs{:});
    else
        [~, edges] = histcounts([data1, data2],name_value_pairs{:});
    end
    % Plot stacked histogram
    [N1, edges1] = histcounts(data1,edges,name_value_pairs{:});
    [N2, edges2] = histcounts(data2,edges,name_value_pairs{:});
    len = length(N1);
    temp = [1:len; 2:len+1]';
    bins = [edges(temp(:,1)); edges(temp(:,2))]';
    bins = mean(bins,2)';
    h1 = mybar([bins; bins]',[N1; N2]','stacked');
    
    
    clear bar_orientation
end



function [h1] = hist_grouped(group,hist_settings)
    %% Plot stacked histogram comparing 2 groups
    
    % Identify all groups with same ctgs value
    ctgs_all = {group.ctgs};
    ctgs_all = cell2mat(ctgs_all);  % Note - this function cannot handle 2D ctgs values at the moment.
    cu = unique(ctgs_all);
    
    legend_arr = get_legendarr(group);
    
    N=length(cu);
    for i = 1:N
        ind = find(ctgs_all == cu(i));
        if N>1; figure; end
        h1 = hist_colors(group(ind),hist_settings);
        legend(legend_arr(ind));
    end
     
end

function [h1] = hist_colors(group,hist_settings)
    
    datastats = {group.datastats};
    datapooled = cell2mat(datastats);
    
    nbins_or_edges = hist_settings.nbins_or_edges;
    name_value_pairs = hist_settings.name_value_pairs;
    
    % Estimate best edges for all data
    if ~isempty(nbins_or_edges); [~, edges] = histcounts(datapooled(:),nbins_or_edges,name_value_pairs{:});
    else [~, edges] = histcounts(datapooled(:),name_value_pairs{:});
    end
    
    % Calculate histogram bins
    N=length(datastats);
    Nbins = zeros(length(edges)-1,N);
    for i = 1:N
        [Nbins(:,i)] = histcounts(datastats{i},edges,name_value_pairs{:});
    end
    
    len = length(edges)-1;
    temp = [1:len; 2:len+1]';
    bins = [edges(temp(:,1)); edges(temp(:,2))]';
    bins = mean(bins,2)';
    bins = repmat(bins(:),1,N);
    h1 = mybar(bins,Nbins,'stacked');
    
end


function [h1] = hist_all(group,plot_type)

    % Nans are due to their being only 1 data point. Remove these by
    % setting to zero (optional)
    
    %nbins_or_edges = hist_settings.nbins_or_edges;
    %name_value_pairs = hist_settings.name_value_pairs;
    %warning('This stuff not yet implemented. Possibly convert hist over to histcounts or just ignore this!');
    
    remove_nans_in_stats = 1;
    
    Ngroups = length(group);
    Nspacings = 100;
    
    pls_all = horzcat(group.data);
    plss_all = horzcat(group.datastats);
    mmin = min(plss_all(:));
    mmax = max(plss_all(:));
    
    plss_spc = linspace(mmin,mmax,Nspacings);
    
    group_n = zeros(Nspacings,Ngroups);
    for i = 1:Ngroups
        [group_n(:,i), xc] = hist(group(i).datastats,plss_spc);    % Old code using hist
    end
    
    
%     group_nc = (group_n);
%     group_nc = cumsum(group_nc,1);
%     group_ncc = cumsum(group_nc,2);    
    
    group_ncc = cumsum(group_n,1);    

    cells_all = horzcat(group.cells);
    
    
    
    switch plot_type
        case 'hist'
            N=length(group);
            for i = 1:N
                subplotsq(N,i);
                h1 = mybar(xc,group_n(:,i));
            end
        case 'histall_imagesc'
            h1 = imagesc(group_n);
            legend_arr = get_legendarr(group,3,0);
            set(gca,'XTickLabel',legend_arr,'XTickLabelRotation',90,'TickLabelInterpreter',group(1).interpreter);
        case 'histall_pdf'
            h1 = plot(xc,group_n);
            %h1 = mybar(xc,group_n);
            legend_arr = get_legendarr(group);
            legend(legend_arr);
        case 'histall_cdf'
            h1 = plot(xc,group_ncc);
            legend_arr = get_legendarr(group);
            legend(legend_arr);
    end
    
    
    
    
    
    
end


function varargout = mybar(varargin)
    % Plot either horizontal or vertical bar graph depending on value of
    % this global variable
    global bar_orientation
    
    if strcmp(bar_orientation,'horiz');
        [varargout{1:nargout}] = barh(varargin{:});
    else
        [varargout{1:nargout}] = bar(varargin{:});
    end
end

