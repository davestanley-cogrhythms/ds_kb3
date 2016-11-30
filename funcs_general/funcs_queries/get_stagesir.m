


function [boundsir, name, nameshort] = get_stagesir(stage, avoid_visual_response, avoid_motor_response)
    % Pull out boundaries of a given stage
    % Boundary outputs are in indices!
    % Boundary outputs are relative!
    % Note that, by convention, all data in a given
    % stage should be >= boundsir(1) and <= boundsir(2).
    
    % Meaning of stages:
    % stage = 1 - pre sample on
    % stage = 2 - sample stage
    % stage = 3 - delay stage
    % stage = 4 - sample + delay stage
    % stage = 5 - all
    % stage = 6 - Test1 - match or non-match
    % stage = 7 - TDelay - test delay
    % stage = 8 - Test2 - match (if prior was non-match)
    
    % Other stages (not yet used)
    % Stage = 0 - fixation
    % stage = -1 - cue
    
    name = 'Undefined';
    nameshort = 'U';
    
    if nargin < 2
        avoid_visual_response = 1; % Jefferson included a 100ms delay to avoid the visual response
    end
    
    if nargin < 3
        avoid_motor_response = 0;
    end

    switch stage
        case -2
            boundsir = [-2.0 -1.0];
            name = 'PreCue';
            nameshort = 'PC';
        case -1
            boundsir = [-1 -0.5];
            name = 'Cue';
            nameshort = 'C';
            
        case 0 
            boundsir = [-0.5 0];
            name = 'Fixation';
            nameshort = 'F';
            
        case 1
            boundsir = [-1.0 0];
            name = 'Pre Sample On';
            
        case 2
            % Sample stage
            boundsir = [0 0.6];
            if avoid_visual_response
                boundsir = [0.1 0.6];    % The value it should be according to Roy et al (i.e. only shift starting point)
                %boundsir = [0.1 0.7];   % The value I was using to get better results (shift starting and ending points)
                                         % Now I've gone back to Roy's value, above
            end
            name = 'Sample';
            nameshort = 'S';
            
        case 3
            % Delay stage
            boundsir = [0.6 1.6];
            %boundsir = [0.9 1.6];
            if avoid_visual_response
                boundsir = [0.9 1.7];   % The values from Roy et al
            end
            name = 'Delay';
            nameshort = 'D';
            
        case 4
            boundsir = [-2.0 1.6];
            if avoid_visual_response
                %boundsir = [-2.0 1.7];
                boundsir = [-1.5 2.0];
            end
            name = 'Pre Sample Delay';
            nameshort = 'A';
            
        case 5
            boundsir = [-2 3.8];
            name = 'All';
            nameshort = 'A';
        case 6
            boundsir = [1.6 2.2];
            if avoid_visual_response
                boundsir = [1.7 2.3];
            end
            if avoid_motor_response
                boundsir = [1.7 1.773];
            end
            name = 'Test';
            nameshort = 'T';
        case 7
            boundsir = [2.2 3.2];
            if avoid_visual_response
                boundsir = [2.5 3.3];
            end
            name = 'TDelay';
            nameshort = 'TD';
        case 8
            boundsir = [3.2 3.8];
            if avoid_visual_response
                boundsir = [3.3 3.9];
            end
            name = 'Test2';
            nameshort = 'T2';
            
    end
    
    dt = get_dt;
    
    boundsir = round(boundsir/dt);
    boundsir(2) = boundsir(2)-1;

end