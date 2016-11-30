


function plot_sfc_vs_time(currspike,curr_stage,ctgsetli,filename,unitname,curr_unit)
   
    
    figure('Position',[240    86   990   712]); 
    unitname_clean = strrep(unitname, '_', ' ');
    %filename = strrep(filename, '.mat', '');
    
    
    %ctgsetli = ctgsetli(:,[2 1 4 3 5:9]);       % REVERSE CTG1 AND CTG2 TO
    %MATCH JEFFERSON! NEVERMIND - I just do it explicitly in the code
    %below.
    
    
    sfc_mode = 2.2;

    [mode_group, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix] = build_sfcmode(sfc_mode, mode_subgroups);
    sfc_mode_group = floor(sfc_mode); % Mode group might be more than 1 dimension, so use floor instead (i.e. sfc_mode = 23.2 would yield mode_group == [2 3] instad of 23)

    sfcfile_sample = fullfile(getpath('path_buffer_curr'),['mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(2)],filename);
    sfcfile_delay = fullfile(getpath('path_buffer_curr'),['mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(3)],filename);


    neither_exist = 1;

    if exist(sfcfile_delay,'file')
        temp = load(sfcfile_delay);
        [h1, h3] = plot_sfc_supporting(temp.sfc,curr_unit,[323 324],'Delay'); clear temp;
        neither_exist = 0;
    end

    if exist(sfcfile_sample,'file')
        temp = load(sfcfile_sample);
        [h1, h3] = plot_sfc_supporting(temp.sfc,curr_unit,[321 322],'Sample'); clear temp;
        neither_exist = 0;
    end


    if neither_exist; plot_sfc_on = 0; end      % Revert to plotting the rasters.

    
    
    % Get ttest results if available:
    schemespath = fullfile(getpath('path_buffer_curr'),'preferred_units');
    schemesname = 'preferred_units_dave.mat';
    
    pschemes = lookup_unit_preferred(filename,unitname_clean,schemespath,schemesname);
    
    
    % Handle leftmost plots - for Ctg1 vs Ctg2
    %subplotset = [.05 .45 .4 .5; .05 .1 .4 .25; ];      % Row 1: Top left ; Row 2: Bottom left
    
    currspike_ctg1 = currspike(:,ctgsetli(:,1)); 
    currspike_ctg2 = currspike(:,ctgsetli(:,2));
    
    if ~plot_sfc_on; h1 = subplot(3,2,[1 3]); plot_raster_conmpare(currspike_ctg1,currspike_ctg2,curr_stage,'r','b');end
    h2 = subplot(3,2,[5]);
    color='r'; [line1] = plot_raster_mean(currspike_ctg1,curr_stage,color);
    color='b'; [line2] = plot_raster_mean(currspike_ctg2,curr_stage,color);
    add_stages; xlabel('Time(ms)');ylabel('Spikes/s');
    legend([line1 line2], 'Ctg1','Ctg2');
    
    
    % Apply title
    titlestr = [filename ' unit ' unitname_clean];
    if ~isempty(pschemes)
        if pschemes(1)
            titlestr = [titlestr ' SampleA'];
        end
        if pschemes(3)
            titlestr = [titlestr ' DelayA'];
        end
    end
    subplot(h1); title(titlestr)
    
    
    % Handle rightmost plots - for Ctg3 vs Ctg4
    %subplotset = [.55 .45 .4 .5; .55 .1 .4 .25; ];      % Row 1: Top right ; Row 2: Bottom right
    currspike_ctg1 = currspike(:,ctgsetli(:,3)); 
    currspike_ctg2 = currspike(:,ctgsetli(:,4));
    
    if ~plot_sfc_on; h3 = subplot(3,2,[2 4]);plot_raster_conmpare(currspike_ctg1,currspike_ctg2,curr_stage,'g','k');end
    h4 = subplot(3,2,[6]);
    color='g'; [line1] = plot_raster_mean(currspike_ctg1,curr_stage,color);
    color='k'; [line2] = plot_raster_mean(currspike_ctg2,curr_stage,color);
    add_stages; xlabel('Time(ms)');ylabel('Spikes/s');
    legend([line1 line2], 'Ctg1','Ctg2');
    
    % Apply title
    titlestr = [];
    if ~isempty(pschemes)
        if pschemes(2)
            titlestr = [titlestr ' SampleB'];
        end
        if pschemes(4)
            titlestr = [titlestr ' DelayB'];
        end
    end
    subplot(h3); title(titlestr)
    
end


function plot_raster_conmpare(currspike_ctg1, currspike_ctg2, curr_stage, colour1, colour2)

    start_row = size(currspike_ctg2,2)+1;
    plot_raster_component(currspike_ctg1,curr_stage,start_row,colour1);

    start_row = 1; 
    plot_raster_component(currspike_ctg2,curr_stage,start_row,colour2);
    
    ticks = [0 size(currspike_ctg2,2) (size(currspike_ctg2,2) + size(currspike_ctg2,1))];
    tick_labels = [0 size(currspike_ctg2,2) size(currspike_ctg1,2)];
    add_stages;
    
    set(gca,'YTick',ticks); set(gca,'YTickLabel',tick_labels);

end


function plot_raster_component(currspike_ctg,curr_stage,start_row,color)

    plot_debug = 0;

    dt = get_dt;
    dt = dt * 1000;
    N = size(currspike_ctg,1);
    %t = (1:N)*dt;
    
    [boundsir] = get_stagesir(curr_stage);
    t = dt * (boundsir(1):boundsir(2));
    
    CS = currspike_ctg';        % Currspike transposed for rasterized plotting
    Ntraces = size(CS,1);

    rasterX = repmat(t,Ntraces,1);
    rasterY = repmat([1:Ntraces]',1,N) + start_row - 1;

    index = find(CS == 1);
    rasterX = rasterX(index); rasterY = rasterY(index);

    if plot_debug
        imagesc([1 N],[1 Ntraces],CS)
    end
    hold on; plot(rasterX,rasterY,[color '.'],'MarkerSize',10);

end


function [line] = plot_raster_mean(currspike_ctg,curr_stage,color)
    
    plot_raw = 0;

    dt = get_dt;
    dt = dt * 1000;
    N = size(currspike_ctg,1);
    %t = (1:N)*dt;
    
    [boundsir] = get_stagesir(curr_stage);
    t = dt * (boundsir(1):boundsir(2));
    
    
    CS = currspike_ctg';        % Currspike transposed for rasterized plotting

    muCS = mean(CS);
    muCS = muCS / get_dt;
    if plot_raw; plot(t,muCS,[color ':']);end
    
    %smCS = sgolayfilt(muCS,3,51);   % Apply 3rd-order filter
    smCS = sgolayfilt(muCS,3,round(151/dt));   % Apply 3rd-order filter, 151 ms 
    %clf;
    hold on; line = plot(t,smCS,color,'LineWidth',1);
    
    ylim([0 max(smCS)]*2);

end






function add_stages

    stages_of_interest = [-1 0 2 3];
    yl = get(gca,'YLim');
    xl = get(gca,'XLim');
    colour = 'kbgrmc';
    for i = 1:length(stages_of_interest)
        stage = stages_of_interest(i);
        [boundsir, name] = get_stagesir(stage);
        %boundsir = boundsir*1000;   % Convert to ms
        
        hold on; plot([boundsir(1) boundsir(1); boundsir(2) boundsir(2)]',[yl(1) yl(2);yl(1) yl(2)]',colour(i),'LineWidth',2);
        text(mean(boundsir),yl(2),name(1),'FontSize',16,'Color',colour(i));
    end
    
    xlim(([xl]-mean(xl))*1.2 + mean(xl));
end



function [h1, h2] = plot_sfc_supporting(sfc,curr_unit,subplots,label)
    
    colours = 'rbgk';
    sfc_ctg1 = sfc{1}.Cave(:,curr_unit);
    sfc_ctg2 = sfc{2}.Cave(:,curr_unit);
    sfc_ctg3 = sfc{3}.Cave(:,curr_unit);
    sfc_ctg4 = sfc{4}.Cave(:,curr_unit);
    f = sfc{1}.f;

    h1 = subplot(subplots(1));
    hold on; plot(f,sfc_ctg1,colours(1));
    hold on; plot(f,sfc_ctg2,colours(2));
    ylabel('SFC'); xlabel('f (Hz)')
    legend(label)
    h2 = subplot(subplots(2));
    hold on; plot(f,sfc_ctg3,colours(3));
    hold on; plot(f,sfc_ctg4,colours(4));
    xlabel('f (Hz)')
    legend(label)

end