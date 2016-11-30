

function fontsz = get_fontlevel(level,scaling_factor)
    % Default font levels for images

    if nargin < 2
        scaling_factor = 1.0;
    end
    
    if nargin < 1
        level = 3;
    end
    
%     if ~is_newgraphics_matlab           % Server Linux running old matlab
%         matlab_version_scale = 1.0;    % Default
%         %matlab_version_scale = 2.5;    % For saving spectrograms
%     end
%     
%     if is_newgraphics_matlab && ~ismac  % Office Linux box
%         matlab_version_scale = 1.0;
%     end
%     
%     if is_newgraphics_matlab && ismac   % Macbook pro
%         matlab_version_scale = 1.2;
%     end
    
    ScreenSize = get(0,'screensize');
    screen_scale = ScreenSize(4)/900;       % Mac screen is 1400x900, so this value by default is 1.0

    switch level
        case 1
            fontsz = 40*scaling_factor*screen_scale;
        case 2
            fontsz = 27*scaling_factor*screen_scale; % Slightly enlarged
        case 3              
            fontsz = 10*scaling_factor*screen_scale; % Matlab default
        case 4              
            fontsz = 20*scaling_factor*screen_scale; % Fixed
    end
    
end