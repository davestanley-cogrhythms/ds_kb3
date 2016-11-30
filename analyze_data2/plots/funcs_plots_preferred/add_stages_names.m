

function add_stages_names (mystages,mycolor)

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
    positions = mean(positions,2);
    
    

    cy = ylim(gca);
    ypos = cy(1) + diff(cy)*9/10;
    
    
    for i = 1:length(positions)
        %text(positions(i),ypos,shortnames{i},'Color',mycolor,'HorizontalAlignment','center','FontSize',26);
        text(positions(i),ypos,shortnames{i},'Color',mycolor,'HorizontalAlignment','center');
        
    end

end