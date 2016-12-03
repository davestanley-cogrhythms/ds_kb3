

%% Main function
% [hsp, out] = plot_spectrogram_custstruct(group,opts)
function [hsp, out] = plot_spectrogram_custstruct(group,opts,overlay_opts,stats_opts)

    % Set defaults
    if nargin < 2
        opts = struct;
    end
    
    if nargin < 3
        overlay_opts = struct;
    end
    
    if nargin < 4
        stats_opts = struct;
    end
    
    abscissa = group.xdata;
    
    % Opts structure
    if ~isfield(opts,'paperfig_mode'); opts.paperfig_mode = 0; end % paperfig_mode: 1=clean plots for paper; 0=plots with additional information for viewing
    if ~isfield(opts,'do_sgolay'); opts.do_sgolay = 0; end
    if ~isfield(opts,'stats_mode'); opts.stats_mode = 0; end
    if ~isfield(opts,'symmetric_axes'); opts.symmetric_axes= 0; end
    if ~isfield(opts,'uniform_zaxis_allgroups'); opts.uniform_zaxis_allgroups = 0; end
    if ~isfield(opts,'do_subplots'); opts.do_subplots = 0; end
    if ~isfield(opts,'max_subplots_per_fig'); opts.max_subplots_per_fig = 16; end
    if ~isfield(opts,'show_range_stats'); opts.show_range_stats = 1; end    % Rectangle associated with stats data
    if ~isfield(opts,'show_range_perm'); opts.show_range_perm = 1; end      % Rectangle associated with permutation data (used to select significant cells)
    if ~isfield(opts,'show_text_stats'); opts.show_text_stats = 1; end
    if ~isfield(opts,'show_text_perm'); opts.show_text_perm = 1; end
    
    % Import variables from opts structure
    paperfig_mode = opts.paperfig_mode;
    do_sgolay = opts.do_sgolay;
    symmetric_axes = opts.symmetric_axes;
    uniform_zaxis_allgroups = opts.uniform_zaxis_allgroups;
    do_subplots = opts.do_subplots;
        max_subplots_per_fig = opts.max_subplots_per_fig;
        use_subplot_grid = 1;
    show_range_stats = opts.show_range_stats;
    show_range_perm = opts.show_range_perm;
    show_text_stats = opts.show_text_stats;
    show_text_perm = opts.show_text_perm;
    
    % Overlay options structure
        % For transparency
    overlay_opts = struct_addDef(overlay_opts,'do_transparency',0);                 % Flag for doing transparency
        % For contours
    overlay_opts = struct_addDef(overlay_opts,'do_contours',0);                     % Flag for doing contours
    overlay_opts = struct_addDef(overlay_opts,'contour_nv',10);    % Contour param n or v - number of contour lines or levels
    overlay_opts = struct_addDef(overlay_opts,'contour_linespec',{'k.'});           % Linespec
    
    % Stats options structure (using new struct_addDef command here)
    stats_opts = struct_addDef(stats_opts,'stats_displaymode',0);                   % Calculates stats on the data and copies it over to contours or overlay or both (data_overlay1 and data_overlay2, respectively).
                                        % 0 = no stats
                                        % 1 = display as transparency
                                        % 2 = display as contours
                                        % 3 = both transparency and contours
    stats_opts = struct_addDef(stats_opts,'statsfunc',@signrank);                          % Stats function to use
    stats_opts = struct_addDef(stats_opts,'stats_comparison',0);                           % Value against which data is compared
    stats_opts = struct_addDef(stats_opts,'contours_alphas',[0.01 0.001 0.0001]);          % Chosen alpha value
    stats_opts = struct_addDef(stats_opts,'transparency_alpha',0.01);                      % Chosen alpha value
    
    
                                    
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
            group(i).data = group(i).data(2:end-1,:,:);
            group(i).xdata = group(i).xdata(2:end-1);
            abscissa = abscissa(2:end-1);
        end
        
    end
    
    % Get max and min z-axis
    if uniform_zaxis_allgroups
        currmax = -Inf;
        currmin = Inf;
        for i = 1:length(group);
            [data] = extract_data(group(i),abscissa,do_sgolay,group(1).xlims_desired,group(1).ylims_desired);
            dat = squeeze(mean(data,2));
            currmin = min(currmin, min(dat(:)));
            currmax = max(currmax, max(dat(:)));
        end
    end

    
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
        
        % Pull out data
        [data, x, x2] = extract_data(group(i),abscissa,do_sgolay,group(1).xlims_desired,group(1).ylims_desired);
        %x2 = group(i).xdata2;
        
        
        % Calculate stats if necessary
        if stats_opts.stats_displaymode > 0
            [group(i).data_overlay1, group(i).data_overlay2, overlay_opts] = calc_pvals(group(i).data,stats_opts,overlay_opts); end
        
        % Plot image
        h=imagesc(x2,x,squeeze(mean(data,2))); set(gca,'YDir','normal');

        % Add overlays - transparency and contours as needed
        draw_contours(group(i),h,group(1).xlims_desired,group(1).ylims_desired,overlay_opts);
        
        legend_arr = get_legendarr(group);

        ylabel(group(1).xdata_name(1,:));
        xlabel(group(1).xdata2_name(1,:));

        % Adjust X and Y axes   
        % This code is no longer needed because I now adjust x/y lims by
        % simply cropping the original data matrix. This is useful because
        % the autoscaling of clims then does not take into account huge
        % values in testing (T) stage.
%         if ~isempty(group(1).xlims_desired); xlim(group(1).xlims_desired); end
%         if ~isempty(group(1).ylims_desired); ylim(group(1).ylims_desired); end
        
        % Adjust Color scale axis
        if uniform_zaxis_allgroups
            caxis([currmin currmax]);
        end
        
        if symmetric_axes       % Make color bar be centered at zero.
            cl = caxis;
            cl = max(abs(cl));
            caxis([-cl cl]);
        end
        if ~isempty(group(i).zlims_desired); caxis(group(i).zlims_desired); end
        
        % Add stage separators
        add_stages_separators([],'w'); add_stages_names([],'w');
        title(legend_arr{i});       % Always show titles
        
        % Set colormap
        if ~isempty(group(i).colourmap_desired);
            colormap(group(i).colourmap_desired);
        else
            colormap jet
        end
        
        % Add colorbar
        if do_subplots;
            %freezeColors
            %cbfreeze(colorbar);
            %hsp.colorbar;
            if uniform_zaxis_allgroups
                % If all colorbars are uniform, only show it for the last
                % subplot
                %if i == length(group); colorbar; end
                colorbar;
            else
                colorbar;
            end
            
            set(gca,'FontSize',10);
        else
            colorbar;
        end
        
        
        
        
        if show_range_perm
            if ~isempty(group(i).full_bandwidth_perm)
                fullband_T = group(i).Nwind_perm*get_dt;
                fullband_F = group(i).full_bandwidth_perm;
                center_T = mean(group(i).timeband_perm);
                center_F = mean(group(i).freqband_perm);
                ctgs = group(i).ctgs;
                mylabel = group(i).tf_label_perm;
                add_spectrogram_tf_rectangle (center_T,center_F,fullband_T,fullband_F,mylabel,'k', '-', show_text_perm,ctgs)
            end
        end
        
        if show_range_stats
            fullband_T = group(i).Nwind*get_dt;
            fullband_F = group(i).full_bandwidth;
            center_T = mean(group(i).timeband_stats);
            center_F = mean(group(i).freqband_stats);
            ctgs = group(i).ctgs;
            mylabel = group(i).tf_label_stats;
            add_spectrogram_tf_rectangle (center_T,center_F,fullband_T,fullband_F,mylabel,'w', '--', show_text_stats,ctgs)
        end
        
        
    end
    
    
    out = [];
end

%% Supporting functions
function Y=ste(varargin)
    Y = std(varargin{:}) / sqrt(size(varargin{1},varargin{3}));
end


function [data_out, xout, x2out] = extract_data(group,abscissa,do_sgolay,xld,yld,myfieldname)

    if nargin < 6
        myfieldname = 'data';
    end
    
    data = group.(myfieldname);
    xdata = group.xdata;
    if ~isempty(xdata); xout = xdata;
    else xout = abscissa;
    end
    
    x2out = group.xdata2;

    % If data still isn't real, we convert it to real, only after
    % taking mean.
    if ~isreal(data)
        data = mean(data,2);
        data = abs(data);
        %data = pls_complex2angle_centpi(data);
        warning('Complex data received! We should be bootstrapping now, meaning that input data should always be real.');
    end
        
    if ~isempty(data) && do_sgolay
        fprintf('Filtering data with SG order 3 filter.\n');
        df = mode(diff(xout));
        sgf = round(20/df);
        if mod(sgf,2) == 0; sgf=sgf+1;end
        data_out = sgolayfilt(data,3,sgf);
    else
        data_out=data;
    end
    
    % Trim trim y-axis (freq) data as needed
    %yld = group.ylims_desired;
    if ~isempty(yld)   % X data holds frequencies, which go along y-axis
        ind = xout >= yld(1) & xout < yld(2);
        xout = xout(ind);
        data_out = data_out(ind,:,:);
    end
    
    % Trim x-axis (time) data as needed
    %xld = group.xlims_desired;
    if ~isempty(xld)            % X-axis data is time, which is stored in xdata2
        ind = x2out > xld(1) & x2out < xld(2);
        x2out = x2out(ind);
        data_out = data_out(:,:,ind);
    end

end


function [amask] = calc_transparency_mask(data,overlay_opts)

    % Pull out structure data
    transparency_mode = overlay_opts.transparency_mode;
    alphamask = overlay_opts.alphamask;
    statsfunc = overlay_opts.statsfunc;
    stats_comparison = overlay_opts.stats_comparison;
    stats_alpha_sig = overlay_opts.stats_alpha_sig;
    
    % Parameters
    test_plot_on = 0;
    
    
    switch transparency_mode
        case 0
            amask = [];
        case 1
            amask = alphamask;
        case {2,3}
            % Calculate stats to get alpha mask
            sz=size(data); datac = mat2cell(data,ones(1,sz(1)),sz(2),ones(1,sz(3)));
%             amask = cellfun(statsfunc,datac);
            statsfunc2 = @(x) statsfunc(x,stats_comparison);
            amask = cellfun(statsfunc2,datac) ;
            amask = squeeze(amask);
            
            if transparency_mode == 2       % Is transparency 
                
            elseif transparency_mode == 3   % Is contour plot
                C = contourc(amask,[stats_pval,stats_pval]);
            end
            
            % Test plotting
            if test_plot_on
                
                
                figure; imagesc(amask); colorbar
                hold on; plot(C(1,:),C(2,:),'w.');
            end
            
            
    end
    
end


function [alpha_mat, co_mat, overlay_opts] = calc_pvals(data,stats_opts,overlay_opts)

    % Plot switches
    test_plot_on = 0;
    
    % Sim switches
    use_for_loop = 1;       % For loop seems to be a bit faster, surprisingly.
    
    % Defaults
    alpha_mat = [];
    co_mat = [];

    % Import data from options structure
    stats_displaymode = stats_opts.stats_displaymode;
                                % 0 = no stats
                                % 1 = display as transparency
                                % 2 = display as contours
                                % 3 = both transparency and contours
    statsfunc = stats_opts.statsfunc;               % Stats function to use
    stats_comparison = stats_opts.stats_comparison; % Value against which data is compared
    contours_alphas = stats_opts.contours_alphas;                 % Chosen alphas for contours
    transparency_alpha = stats_opts.transparency_alpha;           % Chosen alpha value for transparency
    


    % Calculate stats to get alpha mask
    if ~use_for_loop
        sz=size(data); datac = mat2cell(data,ones(1,sz(1)),sz(2),ones(1,sz(3)));
        statsfunc2 = @(x) statsfunc(x,stats_comparison);
        pvals = cellfun(statsfunc2,datac);
    else
        sz=size(data);
        pvals = zeros([sz(1),1,sz(3)]);
        for i = 1:sz(1);
            for j = 1:sz(3)
                pvals(i,1,j) = statsfunc(data(i,:,j),stats_comparison);
            end
        end
    end
    %pvals = squeeze(pvals);        % Don't squeeze because extract_data expects unsqueezed data (3D matrix)

    if any(stats_displaymode == [1,3])         % Display as transparency
        alpha_mat = double(pvals < transparency_alpha);
        %alpha_mat = (alpha_mat + 1) / 2;       % Transparency matrix ranges between 0.5 and 1.0.

        overlay_opts.do_transparency = 1;
    end

    if any(stats_displaymode == [2,3])         % Display as contours
        co_mat = pvals;
        %co_mat = -1*log10(pvals);              % Transform to log scale - tells # zeros after decimal.
        
        overlay_opts.do_contours = 1;
        overlay_opts.contour_nv = contours_alphas;
    end
    
    % Test plotting
    if test_plot_on
        %figure; imagesc(-1*log10(pvals));
        
        figure; h = imagesc(squeeze(mean(data,2))); colorbar; set(gca,'YDir','normal');
        set(h,'AlphaData',squeeze(alpha_mat));
        %hold on; contour(pvals,[0.01,0.01],'w');
        hold on; contour(squeeze(pvals),[0.01, 0.001],'w');
        %hold on; contour(co_mat,'w');


    end
    
end


function draw_contours(group,h,xlims_desired,ylims_desired,overlay_opts)

    % Add transparency
    if overlay_opts.do_transparency
        [alpha_mat] = extract_data(group,group.xdata,0,xlims_desired,ylims_desired,'data_overlay1');
        alpha_mat = squeeze(mean(alpha_mat,2));
        alpha_mat = (alpha_mat + 1) / 2;       % Transparency matrix should range between 0.5 and 1.0.
        set(h,'AlphaData',alpha_mat);
    end
    
    % Add contours
    if overlay_opts.do_contours
        hold on;
        [co_mat, x, x2] = extract_data(group,group.xdata,0,xlims_desired,ylims_desired,'data_overlay2');
        if ~isempty(co_mat)
            co_mat = squeeze(mean(co_mat,2));
            contour(x2,x,co_mat,overlay_opts.contour_nv,overlay_opts.contour_linespec{:});
        else
            warning('Empty contours matrix');
        end
        %contour(x2,x,co_mat);
    end
    
end
