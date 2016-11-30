

function out = colormap_meld(varargin)
    for i = 1:nargin
        colors = varargin(1:end-1);
        N = varargin{end};
    end
    
    out = [];
    for i = 1:length(colors)-1
        out = [out; colormap_meld_pair(colors{i},colors{i+1},N)];
    end
end


function out = colormap_meld_pair(start_color, stop_color,N)


start_color = start_color(:)';
stop_color = stop_color(:)';

diff = stop_color-start_color;
diff = repmat(diff,N,1);

diff_scale = linspace(0,1,N)';
diff_scale = repmat(diff_scale, 1, size(start_color,2));

out = repmat(start_color,N,1);
out = out + diff_scale .* diff;


end