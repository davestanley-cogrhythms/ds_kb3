

function group = query_colourmap_desired(group,group0,middle_colour)
    % This will map the ctg value to a specific colour. 
    
    if nargin < 3
        middle_colour = [];
    end
    
    allctgs = [group.ctgs];
    
    % We want to estimate what the legends should be. Probably a better way
    % to do this...
    group_temp = group;
    group_temp = group.query_legend(group0);
    alllegends = {group_temp.legend};    % Used to get legends
    clear group_temp
    
    for i = 1:length(group)
        if ~isempty(middle_colour)
            out = legend_to_colour_codes(alllegends{i});
            group(i).colourmap_desired = colormap_meld(get_clist(out(2),'matrix'),[0,0,0],get_clist(out(1),'matrix'),64);
        else
            group(i).colourmap_desired = colormap_meld([0,0,0],get_clist(allctgs(i),'matrix'),64);
        end
    end
    
end


function [out] = legend_to_colour_codes(mylegend)

    mylegend = strrep(mylegend,'Cat','1');  % Blue
    mylegend = strrep(mylegend,'Dog','2');  % Green
    mylegend = strrep(mylegend,'Fat','3');  % Red
    mylegend = strrep(mylegend,'Thin','4'); % Cyan
    mylegend = strrep(mylegend,'Goc','3');
    mylegend = strrep(mylegend,'Tad','4');
    
    mylegend = strrep(mylegend,'SchA','1');
    mylegend = strrep(mylegend,'SchB','2');
    
    % Get rid of extra crap
    mylegend = strrep(mylegend,'vs','');
    mylegend = strrep(mylegend,'vs.','');
    mylegend = strrep(mylegend,'/','');
    mylegend = strrep(mylegend,'irr','');
    mylegend = strrep(mylegend,'60%','');    % For mode 7
    mylegend = strrep(mylegend,'100%','');   % For mode 7

    out = str2num(mylegend);
    

end





