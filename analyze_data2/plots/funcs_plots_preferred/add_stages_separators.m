

function add_stages_separators (mystages,mycolor)

    if nargin < 2
        mycolor = 'k';
    end
    
    if nargin < 1
        mystages = [];
    end
    
    if isempty(mystages)
        mystages = [-1 0 2 3 6 7 8];
    end
    
    [positions, names, shortnames] = query_stages(mystages);
    positions = positions(:)';

    cy = ylim(gca);
    cx = xlim(gca);
    
    
    
    for i = 1:length(positions)
        hold on; plot([positions(i) positions(i)], [cy(1) cy(2)], mycolor, 'LineWidth', 4);
    end
    xlim(cx);

end