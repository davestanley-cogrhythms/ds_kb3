

function scaling_factor = get_scaling_factor
    % Default scaling factor for images

%     
%     if ~is_newgraphics_matlab
%         scaling_factor = 1.0;       % Default
%         %scaling_factor = 0.75;      % For saving spectrograms
%     end
%     
%     if is_newgraphics_matlab
%         scaling_factor = 1.0;
%     end
%     
    
    ScreenSize = get(0,'screensize');
    scaling_factor = ScreenSize(4)/900;       % Mac screen is 1400x900, so this value by default is 1.0
    
end