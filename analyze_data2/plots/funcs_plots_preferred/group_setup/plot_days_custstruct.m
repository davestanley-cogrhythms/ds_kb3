

function [h1, out] = plot_days_custstruct(group,funames1D,opts,datenum_numeric)

    % Set defaults
    if nargin < 2
        opts = [];
    end
    
    if ~exist('datenum_numeric','var'); datenum_numeric = []; end
    
    abscissa = group.xdata;

    if ~isfield(opts,'x_axis_mode'); opts.x_axis_mode= 1; end           % Plotting mode. 1-plot vs file number; 2-plot vs time; 3-plot vs time squished (recording sessions ordered chronologically but with gaps deleted)
    if ~isfield(opts,'as_percent'); opts.as_percent= 1; end
    if ~isfield(opts,'uniform_zaxis_allgroups'); opts.uniform_zaxis_allgroups = 1; end
    if ~isfield(opts,'do_subplots'); opts.do_subplots = 1; end
    if ~isfield(opts,'max_subplots_per_fig'); opts.max_subplots_per_fig = 16; end
    
    
    
    
    % Import variables from structures
    x_axis_mode = opts.x_axis_mode;
    as_percent = opts.as_percent;
    uniform_zaxis_allgroups = opts.uniform_zaxis_allgroups;
    do_subplots = opts.do_subplots;
        max_subplots_per_fig = opts.max_subplots_per_fig;
        use_subplot_grid = 1;
        
        
    % Find indices associated with each day
    ids = cellfun(@fname_to_filedate,funames1D);    % ID's
    idsu = unique(ids);                             % ID's unique
    Ndays=length(idsu);
    indDays = {};
    for i = 1:Ndays
        indDays{i} = [ids == idsu(i)]';
    end
    
    % Total number of electrodes per day
    Ntotals = cellfun(@sum,indDays);
    
    % Get cells associated with each group
    Ngroups = length(group);
    for i = 1:Ngroups
        cells{i} = any(group(i).cells,2);
    end
    func = @(x) any(x.cells,2);
    cells = arrayfunu(func,group);
    cells = [cells{:}];
    
    % Calculate number of cells present on each day
    cells_per_day = zeros(Ndays,Ngroups);
    for j = 1:Ngroups
        for i = 1:Ndays
            cells_per_day(i,j) = sum(cells(indDays{i},j));
        end
    end
    
    if as_percent
        cells_per_day = cells_per_day ./ repmat(Ntotals(:),1,Ngroups) * 100;
    end
    
    switch x_axis_mode
        case 1
            x = 1:Ndays;
        case {2}
            if isempty(datenum_numeric); error('datenum_numeric not specified; cant use x_axis_mode = 2; use 1 instead'); end
%             isL = get_ismonkeyL(funames1D);
%             minL = min(idsu(isL));
%             minO = min(idsu(~isL));
            x = datenum_numeric;
            x = x/365;  % Convert to years
            x = x - 2000; % Subtract 2000
            
        case 3
            if isempty(datenum_numeric); error('datenum_numeric not specified; cant use x_axis_mode = 2; use 1 instead'); end
%             isL = get_ismonkeyL(funames1D);
%             minL = min(idsu(isL));
%             minO = min(idsu(~isL));
            
            if x_axis_mode == 3
                x = 1:Ndays;
                [~,ind] = sort(datenum_numeric);
                cells_per_day=cells_per_day(ind,:);
            end
            
            
    end
    
    
    % Get max and min z-axis
    if uniform_zaxis_allgroups
        currmaxz = -Inf;
        currminz = Inf;
        for i = 1:length(group);
            currminz = min(currminz, min(cells_per_day(:,i)));
            currmaxz = max(currmaxz, max(cells_per_day(:,i)));
        end
    end
    
    % Get max and min x-axis
    currmaxx = -Inf;
    currminx = Inf;
    for i = 1:Ngroups
        currminx = min(currminx, min(x(cells_per_day(:,i) > 1)));   % X vale at min entry that has any data
        currmaxx = max(currmaxx, max(x(cells_per_day(:,i) > 1)));   % X vale at max entry that has any data
    end

    % Query legend
    legend_arr = get_legendarr(group);
    
    % Plot the data
    curr_subplots = 1;
    hsp = [];
    max_subplots_per_fig = min(max_subplots_per_fig,length(group));         % Reduce max subplots if length of group is less
    for i = 1:length(group)
        % Manage figure creation
        if do_subplots;
            [hsp, curr_subplots,returns] = new_subplot(max_subplots_per_fig,curr_subplots,hsp,use_subplot_grid);    % Create a new subplot entry or start a new figure if current fig is full.
            if strcmp(returns,'q') || strcmp(returns,'Q'); break; end
        else
            figure;
        end
        
        h1 = bar(x,cells_per_day(:,i));
        
        switch x_axis_mode
            case 1
                xlabel('File number');
            case 2
                xlabel('Date and time');
            case 3; xlabel('File number sorted');
        end
        
        if as_percent
            ylabel('Percent per day');
        else
            ylabel('Number per day');
        end
        
        legend(legend_arr{i});
        
        % Set ylims
        if uniform_zaxis_allgroups
            ylim([currminz currmaxz]);
        end
        
        xrange = [currminx currmaxx];
        xlim(1.1*(xrange - mean(xrange)) + mean(xrange));
    end
    
    out = cells_per_day;
    xout = x;
end
