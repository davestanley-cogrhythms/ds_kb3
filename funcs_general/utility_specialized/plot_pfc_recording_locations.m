


function varargout = plot_pfc_recording_locations(varargin)
    % Same inputs and outputs as normal plotting command

    new_AP = varargin{1};
    new_ML = varargin{2};

    %figure('Name', 'Recording locations PFC', 'Position', [300 300 600 600],'NumberTitle','off');

    %title([FileName(1:end-4), ' PFC  as of ', date])

    hold on
    rectangle('Position',[0,0,25,25],'Curvature',[1,1])
    rectangle('Position',[3,3,19,19],'Curvature',[1,1], 'LineStyle','--')
    [varargout{1:nargout}] = plot(new_AP, new_ML,varargin{3:end});
    text(25, 12.5, '>', 'FontSize', 30)
%     hold off

    %legend('no units', 'units', 'bone')
    xlabel('Post <-------------------> Ant')
    ylabel('Lat <-------------------> Med')

    axis square
    axh = gca;
    set(axh,'XTick',[1:2:25]);
    set(axh,'YTick',[1:2:25]);

    grid
end
