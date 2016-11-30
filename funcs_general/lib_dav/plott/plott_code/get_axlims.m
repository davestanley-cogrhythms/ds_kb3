
function [out] = get_axlims(h)
%     Pulls output is of the from the axis handle and stores them in [out]. If no handle
%     is supplied, gca is used
%     INPUT
%         axis handle (default = gca)
%     
%     OUTPUT
%         out is of the form (xmin xmax ymin ymax) - can be used as input to the axis command
    
    if nargin < 1
        h=gca;
    end
    X = get(gca,'XLim');
    Y = get(gca,'YLim');
    out = [X(:); Y(:)]';

end