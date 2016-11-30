

function h = imagesq(varargin)
    % Imagesc with squeeze

    h = imagesc(squeeze(varargin{1}),varargin{2:end});

end