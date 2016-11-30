

function h = plotsq(varargin)
    % Plot with squeeze

    h = plot(squeeze(varargin{1}),varargin{2:end});

end