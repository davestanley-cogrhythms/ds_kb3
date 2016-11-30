

function [positions, names, shortnames] = query_stages(s,avoid_visual_response,avoid_motor_response)

    if nargin < 2
        avoid_visual_response = 0;
    end
    
    if nargin < 3
        avoid_motor_response = 0;
    end
    
    % Make things into vectors
    if length(avoid_visual_response) == 1; avoid_visual_response = repmat(avoid_visual_response,1,length(s)); end
    if length(avoid_motor_response) == 1; avoid_motor_response = repmat(avoid_motor_response,1,length(s)); end
    
    [positions, names, shortnames] = arrayfunu(@get_stagesir,s,avoid_visual_response,avoid_motor_response);
     
    positions = cell2mat(positions(:)) * get_dt;
    
end
% 
% function [a b c] = temp(x)
%     [a b c] = get_stagesir(x,0);
% end