

function [hsp] = plot_electrode_locations(group,md,funames1D,opts)

    % Set defaults
    if nargin < 4
        opts = struct;
    end
    
%     if ~isfield(opts,'x_axis_mode'); opts.x_axis_mode= 1; end           % Plotting mode. 1-plot vs file number; 2-plot vs time; 3-plot vs time squished (recording sessions ordered chronologically but with gaps deleted)
%     if ~isfield(opts,'as_percent'); opts.as_percent= 1; end
%     if ~isfield(opts,'uniform_zaxis_allgroups'); opts.uniform_zaxis_allgroups = 1; end
    if ~isfield(opts,'do_subplots'); opts.do_subplots = 1; end
    if ~isfield(opts,'do_separate_figures'); opts.do_separate_figures = 0; end
%     if ~isfield(opts,'max_subplots_per_fig'); opts.max_subplots_per_fig = 16; end
    
    
    % Import variables from structures
%     x_axis_mode = opts.x_axis_mode;
%     as_percent = opts.as_percent;
%     uniform_zaxis_allgroups = opts.uniform_zaxis_allgroups;
    do_subplots = opts.do_subplots;
    do_separate_figures = opts.do_separate_figures;
%         max_subplots_per_fig = opts.max_subplots_per_fig;
        use_subplot_grid = 1;
    noise = 0.25;
    
    
    % Pull out electrode coordinates
    temp = [md.electrode_locations];
    AP = [temp.AP];
    ML = [temp.ML];
    
    % Add noise
    AP = AP + unifrnd(-noise,noise,size(AP));
    ML = ML + unifrnd(-noise,noise,size(ML));
        
    % Get ismonkeyL for electrodes (annoying code I've probably written before)
    recording_names_all = {md.recording_name};
    lfp_names_all = [md.lfp_names];
    Nelects_all = arrayfunu(@(x) length(x.lfp_names), md);
    file_names_all = cellfunu(@(x,y) repmat({x},1,y),recording_names_all,Nelects_all);
    file_names_all = horzcat(file_names_all{:});
    files_and_elects = cellfunu(@(x,y) [x ' ' y],file_names_all,lfp_names_all);
    is_monkeyL_elects = get_ismonkeyL(files_and_elects);
    is_monkeyL_elects = is_monkeyL_elects(:);
    
    % Find whether monkeyL and/or monkeyO are present in cells data
    is_monkeyL = get_ismonkeyL(funames1D);
    is_monkeyL = is_monkeyL(:);
    temp=[group.cells];
    allcells = any(temp,2);
    monkey_flags = [ any(is_monkeyL(allcells) == 1) any(is_monkeyL(allcells) == 0)];    % [MonkeyL is present; MonkeyO is present]
    
    
    % List of colors to use
    clist = get_clist;
    clist(1:3) = clist([1,3,2]);    % Make red come first since its a brighter colour; for following paper colour scheme
%     clist(1:3) = clist([3,1,2]);    % Make red come first since its a brighter colour; for following paper colour scheme
    for i = 1:length(group)
        
        % Create new figure and draw all electrode locations in background
        if i == 1 || do_separate_figures
            figl;
            if do_subplots && sum(monkey_flags) > 1;
                hsp = subplot_grid(1,2);
                hsp.set_gca(1); cells_temp = is_monkeyL_elects; hold on; plot_pfc_recording_locations(AP(cells_temp),ML(cells_temp),['k.']);
                hsp.set_gca(2); cells_temp = ~is_monkeyL_elects; hold on; plot_pfc_recording_locations(AP(cells_temp),ML(cells_temp),['k.']);
            else
                inds = any([is_monkeyL_elects.*monkey_flags(1), ~is_monkeyL_elects.*monkey_flags(2) ],2);
                hsp = plot_pfc_recording_locations(AP(inds),ML(inds),['k.']);
            end
        end
    
        % Highlight electrodes in current group (based on group.cells)
        mysymbol='o'; mycolor = clist(i);
        if do_subplots && sum(monkey_flags) > 1;
            hsp.set_gca(1); cells_mask = is_monkeyL;  hold on; plot_recording_group_single(AP,ML,group(i),mycolor,mysymbol,cells_mask);
            hsp.set_gca(2); cells_mask = ~is_monkeyL; hold on; plot_recording_group_single(AP,ML,group(i),mycolor,mysymbol,cells_mask);
        else
            hsp = plot_recording_group_single(AP,ML,group(i),mycolor,mysymbol,[]);
        end
    end
    

end

function h = plot_recording_group_single(AP,ML,group,mycolor,mysymbol,cells_mask)
    % Group has to be length(group) = 1 (i.e. non array)
    
    if nargin < 6
        cells_mask = [];
    end
    
    cells = group.cells;
    cells = any(cells,2);   % Not sure if this is really necessary!
    
    if isempty(cells_mask)
        cells_mask = true(size(cells));
    end
    
    
    if length(cells) > 2*length(AP)       % Cheap hack to guess whether cells are for FFC or PSD; should really save this in the metadata somewhere
        
        % Find allowed inds
        inds = cells & cells_mask;   % Allowed indices cells & cells_mask
        
        % Rebuild full mypairs matrix
        temp = group.metadata.mypairs;
        mypairs = zeros(length(cells),2);
        mypairs(cells,:) = temp;
        
        mypairs = mypairs(inds,:);  % Select allowed pairs
        
        all_elects = [mypairs(:,1); mypairs(:,2)];
        h = plot_pfc_recording_locations(AP(all_elects),ML(all_elects),[mycolor mysymbol],'LineWidth',2);
        %hold on; h2 = plot(AP(mypairs)',ML(mypairs)',['-' mycolor]);
        alpha_transparency = 0.25;
        hold on; h2 = patchline(AP(mypairs)',ML(mypairs)','edgecolor',mycolor,'edgealpha',alpha_transparency);  %Use this command instead to get partial transparency

    elseif length(cells) == length(AP)      % Cells in group are same as electrode cells
        inds = cells & cells_mask;
        h = plot_pfc_recording_locations(AP(inds),ML(inds),[mycolor mysymbol],'LineWidth',2);
    else
        error('asdfadf');
    end

end

