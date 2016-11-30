
function [gr] = group_coords2to4(gr)
    
    sz = [2,2];
    
    for i = 1:length(gr)
        coords1 = gr(i).coords;                         % Gr cycles through coords_alt first, then coords.
        coords2 = gr(i).coords_alt;
        coords_new_temp = [sub2ind(sz,coords1(1),coords1(2)), sub2ind(sz,coords2(1),coords2(2))];   % This means coords_alt will change along columns, coords along rows
        gr(i).coords = coords_new_temp;
        
    end
    
end

% sub2ind maps the coords to linear indices:
% For example
%     sub2ind([2,2],1,1)
%     sub2ind([2,2],2,1)
%     sub2ind([2,2],1,2)
%     sub2ind([2,2],2,2)
%     
% maps to 1,2,3,4

